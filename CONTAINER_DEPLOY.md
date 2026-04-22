# CASCADE Container Deployment Guide

This guide walks through running the **entire** CASCADE evolution pipeline
from a single Docker host using two images:

| Image | Role | GPU? | Where it runs |
|------|------|------|---------------|
| `cascade-controller` | FastAPI dashboard + Vast.ai orchestration + React SPA | CPU-only | Your laptop or a $5/mo VPS |
| `cascade-orchestrator` | PXDesign + Protenix + evolution loop | **GPU (A100)** | One ephemeral Vast.ai instance per parallel run |

Every parallel run gets its own VPS. When the run completes (or fails), the
controller rsyncs the artifacts back, promotes them into the canonical
`outputs/` layout, and destroys the instance.

---

## 1. Prerequisites

- Docker (with `docker compose`), Git, a shell.
- An SSH keypair registered with your Vast.ai account:
  ```bash
  ssh-keygen -t ed25519 -f ~/.ssh/id_ed25519 -N ""
  # add id_ed25519.pub in Vast.ai → Account → SSH Keys
  ```
- A Vast.ai account + CLI key. Log in locally once so the CLI stores the
  token in `~/.config/vastai`:
  ```bash
  pip install --user vastai
  vastai set api-key <YOUR_VAST_API_KEY>
  ```

---

## 2. Build and publish the orchestrator image

The orchestrator image is heavy (~8 GB) because it bakes both the `cascade`
and `pxdesign` conda environments. The repo ships a GitHub Actions workflow
that builds it for you on every push to `main`:

```
.github/workflows/build-orchestrator.yml
  → ghcr.io/<owner>/cascade-orchestrator:latest
```

To trigger it manually:

```bash
git tag orchestrator-v0.1.0
git push --tags
```

If you fork and want a different registry, set
`CASCADE_ORCH_IMAGE=ghcr.io/<you>/cascade-orchestrator:latest` in `.env`.

To build locally (GPU required):

```bash
docker build -t cascade-orchestrator:local -f orchestrator/Dockerfile .
```

See [`orchestrator/README.md`](orchestrator/README.md) for details.

---

## 3. Configure the controller

From the repo root:

```bash
cp .env.example .env
```

Key variables:

```env
CASCADE_ORCH_IMAGE=ghcr.io/twaite11/cascade-orchestrator:latest
VASTAI_CONFIG_DIR=~/.config/vastai
SSH_KEY_DIR=~/.ssh
VAST_SSH_KEY_PATH=/root/.ssh/id_ed25519
```

The controller reads `~/.config/vastai` (read-only bind mount) to authenticate
with the Vast.ai API, and uses the private key at `/root/.ssh/id_ed25519`
inside the container (bind-mounted from your host `~/.ssh`).

---

## 4. Bring up the controller

```bash
docker compose up -d --build
docker compose logs -f controller
```

Then open:

- Dashboard UI: <http://localhost:8000/dashboard>
- API health:   <http://localhost:8000/health>

The controller persists state in these bind mounts:

| Mount | Purpose |
|-------|---------|
| `./outputs` | All runs, logs, promoted artifacts |
| `./metadata` | SQLite catalog, domain metadata, `runs.db` |
| `./data` | Reference DBs, mined hits, HMMs |
| `./jsons` | Per-baseline Protenix JSON inputs |

---

## 5. Click-through walkthrough

**Step 1 — Pick baselines.**
Open the dashboard → `Baselines` tab. Toggle *Validated only* or *Has Phase 1
structure* to narrow the list. Check the enzymes you want to evolve.

**Step 2 — Launch runs.**
Click `Launch N parallel runs`. In the dialog:

- Pick a GPU offer from the Vast.ai offer picker (cheapest A100 80 GB
  instances are suggested by default).
- Set evolution parameters (`--max-generations`, `--variants-per-gen`,
  `--workers`).
- Optionally add a label like `smoke-test`.
- Click `Launch runs`.

The controller:

1. Inserts a `runs` row with status `PROVISIONING`.
2. Calls `vastai create instance <offer_id> <image> --ssh --direct ...`.
3. Polls `vastai show instance` until SSH is reachable.
4. `rsync`es input JSONs + metadata + validated-ID list to the VPS.
5. Spawns an asyncio task that tails `tail -F /workspace/CASCADE/logs/run-<id>.log`
   over SSH and fans the lines out on an in-process `LogHub` — so multiple
   browser tabs can watch the same stream.

**Step 3 — Monitor.**
Jump to the `Runs` tab. Each active run shows:

- live xterm.js panel streaming from the VPS
- cost meter + elapsed timer
- baseline list chips
- `cancel` button (calls `vastai destroy instance`)

**Step 4 — Watch optimized switches appear.**
When the orchestrator emits an optimized switch, the artifact is written to
`/workspace/CASCADE/outputs/runs/<run-id>/optimized_switches/` on the VPS.
On run completion the controller `rsync`es back, promotes the files into
`outputs/optimized_switches/` (prefixed with `<run_id>__` for collision
safety), and invalidates the dashboard service cache. The new variant shows
up in the right-side **Optimized Switches** sidebar within ~15 s.

**Step 5 — Drill in.**
Click any variant ID in the sidebar or the Variants tab to open the detail
drawer: metrics, lineage, crRNA spacer visualization, side-by-side 3D
compare, and artifact download buttons.

Press `?` for the full keyboard-shortcut reference.

---

## 6. Troubleshooting

**`vastai: command not found` inside the controller**

The `vastai` CLI ships in `requirements-controller.txt`. Rebuild with
`docker compose build --no-cache`.

**SSH connect timeouts**

Vast.ai instances take 30–90 s to boot + expose SSH. The controller's
`wait_for_ssh` retries for up to ~5 min. If it still fails, check that your
public key is registered at vastai.com → Account → SSH Keys.

**Logs stop streaming mid-run**

Controller → VPS SSH can drop on flaky networks. The `/api/runs/{id}/logs`
WebSocket endpoint will emit a `[connection closed]` line; refresh the
page to reopen the subscription. Old lines are still in the in-process ring
buffer and replay on reconnect.

**Leftover Vast.ai instances**

Hit `DELETE /api/runs/{id}` (or the red `cancel` button) to destroy a stuck
instance. If the controller itself crashed, run:

```bash
docker compose exec controller vastai show instances
docker compose exec controller vastai destroy instance <id>
```

**Rebuilding after a code change**

```bash
docker compose up -d --build controller
```

---

## 7. Security notes

- The controller is **single-operator**: no auth. Run it on `localhost` or
  behind a VPN / SSH tunnel.
- `~/.ssh` and `~/.config/vastai` are bind-mounted read-only so the
  container cannot tamper with your keys.
- The orchestrator image runs as root inside its VPS. Treat each Vast.ai
  instance as ephemeral; only rsync the artifact tree back, not arbitrary
  files.
