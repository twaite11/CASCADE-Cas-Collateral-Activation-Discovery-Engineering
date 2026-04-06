#!/usr/bin/env bash
set -euo pipefail

HOST="${DASHBOARD_HOST:-0.0.0.0}"
PORT="${DASHBOARD_PORT:-8000}"
export CASCADE_ROOT="${CASCADE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

echo "Starting CASCADE dashboard API on ${HOST}:${PORT}"
echo "UI:    http://${HOST}:${PORT}/dashboard"
echo "Health: http://${HOST}:${PORT}/health"
echo "CASCADE_ROOT=${CASCADE_ROOT}"

uvicorn dashboard_backend.main:app --host "${HOST}" --port "${PORT}"
