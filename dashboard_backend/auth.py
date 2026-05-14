"""API-key authentication for mutating dashboard endpoints (C-2).

CASCADE's dashboard exposes a handful of endpoints that spin up paid Vast.ai
GPU instances, destroy them, and tail their logs.  Before this module, any
caller reachable on port 8000 (which the legacy Dockerfile binds to
``0.0.0.0`` with ``--forwarded-allow-ips "*"``) could:

  * ``POST /api/runs``           - burn the operator's GPU wallet
  * ``DELETE /api/runs/{id}``    - cancel an active scientific run
  * ``WS   /api/runs/{id}/logs`` - tap live job output
  * ``GET  /api/vast/offers``    - enumerate available GPUs (rate-limit foot-gun)

We gate all four behind an opt-in API-key check.  By default (no key set)
the backend is **closed**: every mutating endpoint returns 401 and the only
way to enable them is to deploy with ``CASCADE_API_KEY`` set.  The legacy
"wide-open" behaviour is still reachable for one-off local development by
setting ``CASCADE_API_KEY=disabled`` (we surface a logged warning every time
this happens so it shows up in operator dashboards).

Read-only endpoints (`/api/variants`, `/api/baselines`, ...) are intentionally
not authenticated -- the React UI fetches them anonymously and they expose no
mutating capability or secrets.
"""
from __future__ import annotations

import logging
import os
import secrets
from typing import Optional

from fastapi import Depends, Header, HTTPException, Query, WebSocket, status

log = logging.getLogger(__name__)

API_KEY_ENV = "CASCADE_API_KEY"
HEADER_NAME = "X-Cascade-Api-Key"
QUERY_PARAM = "api_key"  # WebSockets cannot send custom headers from the browser
DISABLED_TOKEN = "disabled"


def _expected_key() -> Optional[str]:
    """Return the configured key, or None if auth is not configured."""
    val = os.environ.get(API_KEY_ENV, "").strip()
    return val or None


def is_auth_disabled() -> bool:
    """True iff the operator has explicitly opted out of API-key gating."""
    return _expected_key() == DISABLED_TOKEN


def _warn_disabled_once() -> None:
    """Log a (single) warning when the operator runs with auth disabled."""
    if not getattr(_warn_disabled_once, "_warned", False):
        log.warning(
            "%s=%r — mutating dashboard endpoints are UNAUTHENTICATED. "
            "This is acceptable for `cascade dev` on localhost but should "
            "never reach production / public networks.",
            API_KEY_ENV, DISABLED_TOKEN,
        )
        _warn_disabled_once._warned = True  # type: ignore[attr-defined]


def require_api_key(
    x_cascade_api_key: Optional[str] = Header(default=None, alias=HEADER_NAME),
    api_key_query: Optional[str] = Query(default=None, alias=QUERY_PARAM),
) -> str:
    """FastAPI dependency for HTTP routes.

    Accept the key either as the ``X-Cascade-Api-Key`` header (preferred) or
    as the ``?api_key=`` query parameter (for browser GETs and curl one-liners
    that can't easily set custom headers).  Comparison is constant-time.
    """
    expected = _expected_key()
    if expected is None:
        # Fail-closed by default: the operator must set the env var.
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=(
                f"{API_KEY_ENV} not configured; mutating endpoints disabled. "
                f"Set the env var (or =disabled for localhost dev) to enable."
            ),
            headers={"WWW-Authenticate": HEADER_NAME},
        )
    if expected == DISABLED_TOKEN:
        _warn_disabled_once()
        return DISABLED_TOKEN

    presented = x_cascade_api_key or api_key_query or ""
    if not presented or not secrets.compare_digest(presented, expected):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="invalid or missing API key",
            headers={"WWW-Authenticate": HEADER_NAME},
        )
    return presented


async def require_api_key_ws(websocket: WebSocket) -> bool:
    """WebSocket equivalent of ``require_api_key``.

    Returns True on success.  On failure it closes the socket with a 4401
    application-level close code and returns False.  Browsers can't easily
    send custom headers when opening a WS, so we accept the key from either
    the ``X-Cascade-Api-Key`` header or the ``?api_key=`` query string.
    """
    expected = _expected_key()
    if expected is None:
        await websocket.close(code=4401, reason="api key not configured")
        return False
    if expected == DISABLED_TOKEN:
        _warn_disabled_once()
        return True

    presented = (
        websocket.headers.get(HEADER_NAME.lower())
        or websocket.query_params.get(QUERY_PARAM)
        or ""
    )
    if not presented or not secrets.compare_digest(presented, expected):
        await websocket.close(code=4401, reason="invalid api key")
        return False
    return True
