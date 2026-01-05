#!/bin/bash

# ============================================
# Confgen Webapp - Celery Startup Script
# (Apache2 handles the web server)
# ============================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# Default to development, override with: ./start.sh prod
ENV_MODE="${1:-dev}"

echo -e "${CYAN}"
echo "╔═══════════════════════════════════════╗"
echo "║     Confgen - Start Celery Worker     ║"
echo "╚═══════════════════════════════════════╝"
echo -e "${NC}"

# Set Flask environment
if [ "$ENV_MODE" = "prod" ]; then
    export FLASK_ENV=production
    echo -e "  ${YELLOW}Mode: PRODUCTION${NC}"
else
    export FLASK_ENV=development
    echo -e "  ${GREEN}Mode: DEVELOPMENT${NC}"
fi
echo ""

# Configuration
APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$APP_DIR/.venv"
LOG_DIR="$APP_DIR/logs"
PID_DIR="$APP_DIR/pids"

mkdir -p "$LOG_DIR" "$PID_DIR"

# Activate virtual environment
if [ -d "$VENV_DIR" ]; then
    echo -e "${YELLOW}[1/3]${NC} Activating virtual environment..."
    source "$VENV_DIR/bin/activate"
    echo -e "  ${GREEN}✓ Done${NC}"
else
    echo -e "${RED}Error: Virtual environment not found at $VENV_DIR${NC}"
    exit 1
fi

# Check Redis
echo -e "${YELLOW}[2/3]${NC} Checking Redis..."
if redis-cli ping &> /dev/null; then
    echo -e "  ${GREEN}✓ Redis is running${NC}"
else
    echo -e "  ${YELLOW}Starting Redis...${NC}"
    sudo systemctl start redis-server 2>/dev/null || redis-server --daemonize yes
    sleep 1
    if redis-cli ping &> /dev/null; then
        echo -e "  ${GREEN}✓ Redis started${NC}"
    else
        echo -e "  ${RED}✗ Failed to start Redis${NC}"
        exit 1
    fi
fi

# Start Celery worker
echo -e "${YELLOW}[3/3]${NC} Starting Celery worker..."
pkill -f "celery -A app.tasks:celery worker" 2>/dev/null
sleep 1

cd "$APP_DIR"
FLASK_ENV=$FLASK_ENV celery -A app.tasks:celery worker --loglevel=info \
    --logfile="$LOG_DIR/celery.log" \
    --pidfile="$PID_DIR/celery.pid" \
    --detach

sleep 3
if [ -f "$PID_DIR/celery.pid" ]; then
    PID=$(cat "$PID_DIR/celery.pid")
    if kill -0 "$PID" 2>/dev/null; then
        echo -e "  ${GREEN}✓ Celery started (PID: $PID)${NC}"
    else
        echo -e "  ${RED}✗ Failed to start Celery. Check $LOG_DIR/celery.log${NC}"
        exit 1
    fi
else
    # Check if process is running anyway
    if pgrep -f "celery -A app.tasks:celery worker" > /dev/null; then
        echo -e "  ${GREEN}✓ Celery started${NC}"
    else
        echo -e "  ${RED}✗ Failed to start Celery. Check $LOG_DIR/celery.log${NC}"
        exit 1
    fi
fi

echo ""
echo -e "${GREEN}════════════════════════════════════════${NC}"
echo -e "${GREEN}  Celery worker is running!${NC}"
echo -e "${GREEN}════════════════════════════════════════${NC}"
echo ""
echo -e "  ${CYAN}Logs:${NC}   $LOG_DIR/celery.log"
echo -e "  ${CYAN}Stop:${NC}   ./stop.sh"
echo ""
echo -e "  ${YELLOW}Usage:${NC}  ./start.sh        (development)"
echo -e "          ./start.sh prod   (production)"
echo ""
