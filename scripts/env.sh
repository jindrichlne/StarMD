# Source this to keep every tool cache (npm, pnpm, XDG) inside the project folder, under .tools/.
_starmd="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
export npm_config_cache="$_starmd/.tools/npm-cache"
export XDG_CACHE_HOME="$_starmd/.tools/xdg/cache" XDG_DATA_HOME="$_starmd/.tools/xdg/data" XDG_STATE_HOME="$_starmd/.tools/xdg/state" XDG_CONFIG_HOME="$_starmd/.tools/xdg/config"
export PNPM_HOME="$_starmd/.tools/pnpm-home"
export npm_config_store_dir="$_starmd/.tools/pnpm-store"
export PUPPETEER_SKIP_DOWNLOAD=1 PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD=1 ASTRO_TELEMETRY_DISABLED=1 npm_config_update_notifier=false
export PATH="$_starmd/.tools/node_modules/.bin:$PATH"
unset _starmd
