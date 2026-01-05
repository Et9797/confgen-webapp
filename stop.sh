#!/bin/bash

# ============================================
# Confgen Webapp - Stop Celery Script
# ============================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${CYAN}"
echo "╔═══════════════════════════════════════╗"
echo "║     Confgen - Stop Celery Worker      ║"
echo "╚═══════════════════════════════════════╝"
echo -e "${NC}"

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PID_DIR="$APP_DIR/pids"

echo -e "${YELLOW}Stopping Celery worker...${NC}"
if [ -f "$PID_DIR/celery.pid" ]; then
    PID=$(cat "$PID_DIR/celery.pid")
    if kill -0 "$PID" 2>/dev/null; then
        kill "$PID"
        rm -f "$PID_DIR/celery.pid"
        echo -e "  ${GREEN}✓ Celery stopped (PID: $PID)${NC}"
    else
        rm -f "$PID_DIR/celery.pid"
        echo -e "  ${YELLOW}✓ Process was not running${NC}"
    fi
else
    pkill -f "celery -A app.tasks:celery worker" 2>/dev/null
    echo -e "  ${GREEN}✓ Stopped${NC}"
fi

echo ""
echo -e "${GREEN}Celery worker stopped.${NC}"
echo ""
