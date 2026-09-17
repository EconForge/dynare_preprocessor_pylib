#!/usr/bin/env bash
# Execute commands inside Windows guest over SSH
set -e

KEY_FILE="$HOME/.win-dev/id_ed25519"
SSH_PORT="2222"
SSH_USER="Docker"
SSH_HOST="localhost"

SSH_OPTS=(
  -i "$KEY_FILE"
  -p "$SSH_PORT"
  -o StrictHostKeyChecking=no
  -o UserKnownHostsFile=/dev/null
  -o LogLevel=ERROR
  -o ConnectTimeout=5
)

SYNC_CMD='robocopy \\host.lan\Data C:\preprocessor /MIR /XD .git .pixi build build-win win32 test_llvm /NFL /NDL /NJH /NJS /nc /ns /np /R:1 /W:1'

# If no arguments given, open an interactive shell in C:\preprocessor
if [ $# -eq 0 ]; then
  ssh -t "${SSH_OPTS[@]}" "${SSH_USER}@${SSH_HOST}" "cmd.exe /k \"${SYNC_CMD} & cd /d C:\\preprocessor\""
else
  ssh "${SSH_OPTS[@]}" "${SSH_USER}@${SSH_HOST}" "cmd.exe /c \"${SYNC_CMD} & cd /d C:\\preprocessor && $*\""
fi
