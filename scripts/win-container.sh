#!/usr/bin/env bash
set -e

CONTAINER_NAME="win-dev"
STORAGE_DIR="$HOME/.win-dev/storage"
OEM_DIR="$HOME/.win-dev/oem"
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

case "$1" in
  start)
    mkdir -p "$STORAGE_DIR" "$OEM_DIR"
    if docker ps --filter "name=^/${CONTAINER_NAME}$" --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
      echo "Container ${CONTAINER_NAME} is already running."
    elif docker ps -a --filter "name=^/${CONTAINER_NAME}$" --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
      echo "Starting existing container ${CONTAINER_NAME}..."
      docker start "${CONTAINER_NAME}"
    else
      echo "Creating and launching ${CONTAINER_NAME}..."
      docker run -d \
        --name "${CONTAINER_NAME}" \
        --device /dev/kvm \
        --cap-add NET_ADMIN \
        -p 8006:8006 \
        -p 2222:22 \
        -p 3389:3389 \
        -v "${SRC_DIR}:/shared" \
        -v "${STORAGE_DIR}:/storage" \
        -v "${OEM_DIR}:/oem" \
        -e VERSION="tiny11" \
        -e RAM_SIZE="8G" \
        -e CPU_CORES="4" \
        -e DISK_SIZE="64G" \
        dockurr/windows
      echo "Container launched. Windows setup has started."
      echo "View progress at http://localhost:8006 or via 'bash scripts/win-container.sh logs'."
    fi
    ;;
  stop)
    echo "Stopping ${CONTAINER_NAME}..."
    docker stop "${CONTAINER_NAME}"
    ;;
  status)
    docker ps -a --filter "name=^/${CONTAINER_NAME}$"
    ;;
  logs)
    docker logs -f "${CONTAINER_NAME}"
    ;;
  destroy)
    echo "Stopping and removing ${CONTAINER_NAME}..."
    docker rm -f "${CONTAINER_NAME}" || true
    ;;
  *)
    echo "Usage: $0 {start|stop|status|logs|destroy}"
    exit 1
    ;;
esac

