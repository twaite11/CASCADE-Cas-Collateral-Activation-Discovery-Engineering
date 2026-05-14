"""Deterministic seeding helpers for CASCADE pipelines.

B-12 fix: previously the evolution orchestrator passed `seed=42` *only* to
the mismatch sequence generator.  Every other source of randomness
(`numpy.random`, `random`, `torch`, `cudnn`) was implicitly unseeded, so
re-running the same Generation N produced different variants, different
HEPN-shift values, and different elite picks.  This module is the single
place where we seed every RNG that matters.

Usage at the top of any orchestrator-like script:

    from utils.determinism import seed_all, resolve_seed
    seed = resolve_seed()                 # picks up --seed / CASCADE_SEED / 42
    seed_all(seed)                        # applies to random / numpy / torch
"""
from __future__ import annotations

import logging
import os
import random
from typing import Optional

log = logging.getLogger(__name__)

DEFAULT_SEED = 42
SEED_ENV = "CASCADE_SEED"


def resolve_seed(cli_value: Optional[int] = None) -> int:
    """Pick the effective seed.

    Precedence (highest first):
      1. ``cli_value`` -- whatever ``--seed`` passed on the command line
      2. ``CASCADE_SEED`` environment variable
      3. ``DEFAULT_SEED`` (42 -- the canonical "deterministic" sentinel)
    """
    if cli_value is not None:
        return int(cli_value)
    env_val = os.environ.get(SEED_ENV, "").strip()
    if env_val:
        try:
            return int(env_val)
        except ValueError:
            log.warning(
                "%s=%r is not an integer; falling back to default %d",
                SEED_ENV, env_val, DEFAULT_SEED,
            )
    return DEFAULT_SEED


def seed_all(seed: int) -> None:
    """Apply ``seed`` to every RNG that affects CASCADE outputs.

    Always seeds ``random`` and ``numpy``.  Conditionally seeds ``torch``
    (CPU + CUDA) and sets cuDNN to deterministic mode if PyTorch is
    importable.  Safe to call repeatedly.
    """
    random.seed(seed)
    try:
        import numpy as np  # noqa: PLC0415
        np.random.seed(seed)
    except ImportError:
        pass

    try:  # PyTorch is optional at dev time but required on the GPU host.
        import torch  # noqa: PLC0415
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        # cudnn.deterministic trades a small perf hit for bit-exact
        # convolution outputs across runs.  Worth it for an evolution
        # pipeline whose entire correctness story depends on
        # reproducibility.
        try:
            import torch.backends.cudnn as cudnn  # noqa: PLC0415
            cudnn.deterministic = True
            cudnn.benchmark = False
        except (ImportError, AttributeError):
            pass
    except ImportError:
        pass

    # Some downstream tools (e.g. MPNN, RFDiffusion via pxdesign) read
    # PYTHONHASHSEED to deduplicate cached prompts.  Set it for the
    # current process so subprocess children inherit a stable value.
    os.environ.setdefault("PYTHONHASHSEED", str(seed))
    log.info("Deterministic mode: seeded random/numpy/torch with seed=%d", seed)
