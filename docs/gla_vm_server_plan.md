# Glasgow VM Server — Deployment Plan

> **Status**: Deployed (`3d8a62c`, 2026-02-20). This is a record of the VM setup as it was done, so the live server may have drifted from it. The VM uses **micromamba** (`~/micromamba`), and the systemd and SELinux paths below depend on it, so it is kept deliberately. `setup_all.sh` still detects micromamba first. For local installs, use mamba/conda as in the README.

Deploy the smFISHProbeDesign Streamlit app on an Oracle Linux 9.6 university VM.

---

## VM Specifications

| Resource | Value |
|----------|-------|
| OS | Oracle Linux Server 9.6 (RHEL-family, `el9`) |
| CPU | Intel Xeon Silver 4310 @ 2.10 GHz — 2 cores / 4 threads |
| RAM | 3.4 GB total |
| Swap | 3.9 GB |
| Root disk (`/`) | 70 GB (61 GB free) |
| Home disk (`/home`) | 425 GB (421 GB free) |

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Package manager | micromamba | Fast, lightweight, no Anaconda licensing |
| Install location | `~/Github/smFISHProbeDesign` | User's home directory, easy `git pull` updates |
| Process manager | systemd | Auto-restart, survives reboots, journal logging |
| Network exposure | IT's HTTPS Nginx config + ProbeDesign location block | IT pre-installed nginx with Let's Encrypt cert and HTTP→HTTPS redirect; we add our app as a location block inside their HTTPS server |
| RepeatMasker | Skip for now | Users can upload pre-masked files; add later if needed |
| Swap | Keep existing 3.9 GB | Already sufficient — genome masking peaks ~3.5 GB |

## Disk Space Budget

| Component | Size | Location |
|-----------|------|----------|
| Repository (code + pseudogene FASTAs) | ~200 MB | `~/Github/smFISHProbeDesign` |
| Conda environment (Python, Bowtie, Streamlit) | ~500 MB | `~/micromamba/envs/probedesign` |
| Pseudogene bowtie indices (3 species) | ~85 MB | `~/Github/smFISHProbeDesign/bowtie_indexes` |
| Genome indices (human + mouse + drosophila) | ~6.7 GB | `~/Github/smFISHProbeDesign/bowtie_indexes` |
| **Total** | **~7.5 GB** | All under `/home` (421 GB free) |

---

## Step-by-step Deployment

### Step 1 — SSH into the VM

```bash
ssh <your-username>@<vm-ip-address>
```

Verify the OS and resources:

```bash
cat /etc/os-release          # Should show Oracle Linux 9.6
free -h                      # Confirm ~3.4 GB RAM + ~3.9 GB swap
df -h /home                  # Confirm ~421 GB free
```

---

### Step 2 — Install System Prerequisites

Install compilers and tools needed by conda packages and bowtie index building.

```bash
sudo dnf install -y git tar unzip wget curl gcc gcc-c++ make bzip2
```

Verify git:

```bash
git --version
```

---

### Step 3 — Install micromamba

micromamba is a fast, standalone conda package manager. Install it into your home directory (no root required for this step).

```bash
# Download and install micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

During install, accept the defaults:
- Install location: `~/micromamba` (or `~/.local/bin/micromamba`)
- Shell init: Yes (adds to `~/.bashrc`)

After install, reload your shell:

```bash
source ~/.bashrc
micromamba --version
```

Expected output: `2.x.x` or similar.

**Configure channels** (one-time):

```bash
micromamba config append channels conda-forge
micromamba config append channels bioconda
micromamba config set channel_priority strict
```

---

### Step 4 — Clone the Repository

```bash
mkdir -p ~/Github
cd ~/Github
git clone https://github.com/JeffLeeLab/smFISHProbeDesign.git
cd smFISHProbeDesign
```

---

### Step 5 — Run the Setup Script Inside a screen Session

The setup script downloads ~6.7 GB of genome indices and takes 15–30+ minutes depending on network speed. **Run it inside `screen`** so the process survives SSH disconnection. If your connection drops mid-download, you can reconnect and reattach to find it still running.

**5a. Install `screen`** (if not already present):

```bash
sudo dnf install -y screen
```

**5b. Start a named screen session:**

```bash
screen -S probedesign_setup
```

You are now inside a persistent session. The terminal looks the same — the difference is the session lives on the server independently of your SSH connection.

**5c. Run the setup script:**

```bash
cd ~/Github/smFISHProbeDesign
chmod +x setup_all.sh
./setup_all.sh
```

**Essential `screen` commands:**

| Action | Keys |
|--------|------|
| Detach (leave running) | `Ctrl-A` then `D` |
| Reattach after reconnecting | `screen -r probedesign_setup` |
| List all sessions | `screen -ls` |
| Kill session when done | `exit` (inside session) |

If your SSH connection drops, simply SSH back in and run `screen -r probedesign_setup` to reattach.

**What the script does:**

1. **Create `probedesign` conda environment** (~500 MB download)
   - Python 3.11, Bowtie 1.3.1, Click, Streamlit, Pandas
2. **Install the `probedesign` package** (editable mode)
3. **Build pseudogene indices** from shipped FASTAs (~1 min)
   - `humanPseudo`, `mousePseudo`, `drosophilaPseudo`
4. **Download genome indices** from AWS (~6.7 GB, ~12–30 min depending on network)
   - Human GRCh38 (~3.5 GB), Mouse mm10 (~3.1 GB), Drosophila BDGP6 (~176 MB)
5. **Run validation tests** (7 checks)

**If you only need pseudogene masking** (skip the 6.7 GB genome download):

```bash
./setup_all.sh --skip-genome
```

You can always re-run `./setup_all.sh` later — it skips already-completed steps, so it is safe to resume an interrupted run.

**Expected output** (all tests pass):

```
Tests: 7/7 passed
All checks passed!
```

When done, exit the screen session:

```bash
exit
```

---

### Step 6 — Verify the Installation

Activate the environment and run quick checks.

```bash
micromamba activate probedesign

# Check tools are available
python --version           # Python 3.11.x
bowtie --version           # bowtie-align-s version 1.3.1
streamlit --version        # 1.x.x
probedesign --version      # 0.1.0
```

**Verify bowtie index auto-detection:**

Both the CLI and the Streamlit app derive the index directory automatically from the script's own location using `Path(__file__).resolve()`, not from the current working directory. This means the path is always `~/Github/smFISHProbeDesign/bowtie_indexes` regardless of where you `cd` before running.

Confirm the indices are where the app expects them:

```bash
ls ~/Github/smFISHProbeDesign/bowtie_indexes/*.1.ebwt 2>/dev/null | head -5
# Should list: humanPseudo.1.ebwt, mousePseudo.1.ebwt, drosophilaPseudo.1.ebwt

ls ~/Github/smFISHProbeDesign/bowtie_indexes/*.1.bt2 2>/dev/null | head -5
# Should list: GCA_000001405.15_GRCh38_no_alt_analysis_set.1.bt2, mm10.1.bt2, drosophila.1.bt2
```

Test CLI probe design with pseudogene masking to confirm the index is found:

```bash
cd ~/Github/smFISHProbeDesign

# No masking — fast sanity check
probedesign design tests/regression/CDKN1A_32/CDKN1A.fa --probes 32 \
  --repeatmask-file tests/regression/CDKN1A_32/CDKN1A_repeatmasked.fa \
  -o /tmp/test_cdkn1a --quiet
cat /tmp/test_cdkn1a_oligos.txt | wc -l   # Should print: 32

# With pseudogene masking (verifies index auto-detection)
probedesign design tests/regression/KRT19_withUTRs/KRT19_withUTRs.fa --probes 32 \
  --pseudogene-mask --quiet -o /tmp/test_krt19
echo "Exit code: $?"   # Should print: Exit code: 0
```

In the Streamlit UI, the "Bowtie index directory" sidebar field is pre-filled with the same auto-detected path. The `check_prerequisites()` function validates the index files exist before any design run and shows a warning banner if anything is missing — it will never silently fail.

**Quick Streamlit smoke test** (verify it starts, then Ctrl+C):

```bash
streamlit run streamlit_app/app.py --server.port 8501 --server.address 0.0.0.0
# You should see: "You can now view your Streamlit app in your browser."
# Press Ctrl+C to stop — the systemd service will manage it from here
```

---

### Step 7 — Create a systemd Service

A systemd unit ensures the app starts on boot, restarts on crash, and logs to journald.

**7a. Find the full path to Streamlit in the conda environment:**

```bash
# While probedesign env is activated:
which streamlit
# Example output: /home/jeff/micromamba/envs/probedesign/bin/streamlit
```

Note this path — you will use it in the service file below.

**7b. Fix SELinux to allow execution from the home directory:**

Oracle Linux 9 (RHEL family) runs SELinux in enforcing mode. systemd (`init_t`) cannot execute binaries under `/home` unless they carry the `bin_t` label. Run this once:

```bash
sudo dnf install -y policycoreutils-python-utils   # provides semanage

sudo semanage fcontext -a -t bin_t "/home/jeff/micromamba(/.*)?"
sudo restorecon -Rv /home/jeff/micromamba/
```

Verify the label was applied:

```bash
ls -Z /home/jeff/micromamba/envs/probedesign/bin/streamlit
# Should show: system_u:object_r:bin_t:s0
```

**7c. Create the service file:**

> **Note — Linux paths are case-sensitive.** Verify your clone directory before pasting: `ls ~/Git*/`

```bash
sudo tee /etc/systemd/system/probedesign.service > /dev/null << 'SERVICEEOF'
[Unit]
Description=ProbeDesign Streamlit Web App
After=network.target

[Service]
Type=simple
User=jeff
Group=jeff

# Use absolute path to app.py — avoids SELinux WorkingDirectory issues
ExecStart=/home/jeff/micromamba/envs/probedesign/bin/streamlit run \
    /home/jeff/Github/smFISHProbeDesign/streamlit_app/app.py \
    --server.port 8501 \
    --server.address 127.0.0.1 \
    --server.baseUrlPath /probedesign \
    --server.headless true \
    --browser.gatherUsageStats false

# Restart policy
Restart=on-failure
RestartSec=5

# Environment — ensure the conda env's bin is in PATH so bowtie is found
Environment="PATH=/home/jeff/micromamba/envs/probedesign/bin:/usr/local/bin:/usr/bin:/bin"
Environment="HOME=/home/jeff"

# Resource limits
LimitNOFILE=65536

# Logging
StandardOutput=journal
StandardError=journal
SyslogIdentifier=probedesign

[Install]
WantedBy=multi-user.target
SERVICEEOF
```

**7d. Enable and start the service:**

```bash
sudo systemctl daemon-reload
sudo systemctl enable probedesign        # Start on boot
sudo systemctl start probedesign         # Start now
```

**7e. Verify the service is running:**

```bash
sudo systemctl status probedesign
```

Expected output should show `active (running)`. If it fails, check logs:

```bash
sudo journalctl -u probedesign -f                   # Follow live logs
sudo journalctl -u probedesign --no-pager -n 50     # Last 50 lines
```

---

### Step 8 — Configure Nginx

Nginx acts as a reverse proxy, forwarding browser requests to the Streamlit process on `127.0.0.1:8501`.

> **On this VM (and likely on any university-managed VM), IT may have already installed and configured Nginx** — including HTTPS certificates, HTTP→HTTPS redirects, and location blocks for other apps (e.g. R Shiny). You must work with that existing config, not replace it. Steps 8a–8c below show how to check and then integrate.

**8a. Install Nginx (if not already present):**

```bash
sudo dnf install -y nginx
sudo systemctl enable --now nginx
```

If `nginx -v` already returns a version, it's installed — skip this step.

**8b. Check for an existing nginx config from IT:**

```bash
ls /etc/nginx/conf.d/
```

Look at any `.conf` files you didn't create. Inspect them:

```bash
sudo grep -rn "ssl\|https\|return 301\|server_name\|listen" /etc/nginx/conf.d/
```

**What you're looking for:**

| What you see | What it means |
|---|---|
| `return 301 https://` on port 80 | IT is redirecting ALL HTTP to HTTPS — a port-80-only config will never be reached |
| `listen 443 ssl` with a `server_name` matching the VM hostname | IT has an HTTPS server block; this is where your location blocks must go |
| `ssl_certificate /etc/letsencrypt/...` | Let's Encrypt cert already provisioned — HTTPS is fully working |
| `proxy_pass http://127.0.0.1:3838` | IT is running R Shiny Server; your app shares the same Nginx |

**On this VM** (`shinytest.mvls.gla.ac.uk`), IT's file is `/etc/nginx/conf.d/shiny.conf` and it does exactly this: redirects HTTP→HTTPS on port 80, and proxies `location /` to R Shiny on port 3838 over HTTPS with a Let's Encrypt cert.

**8c. Create `streamlit_apps.conf` with only the upstream definition:**

This file defines the named upstream for Streamlit (http-context, not inside a server block). The actual server and location blocks go in IT's HTTPS config (Step 8d).

```bash
sudo tee /etc/nginx/conf.d/streamlit_apps.conf > /dev/null << 'EOF'
# ProbeDesign Streamlit upstream
upstream probedesign {
    server 127.0.0.1:8501;
}

# Future apps — add an upstream per app:
# upstream rnaviewer {
#     server 127.0.0.1:8502;
# }
EOF
```

**8d. Add ProbeDesign location blocks to IT's HTTPS server block:**

Open IT's nginx config (e.g. `shiny.conf`) and add the ProbeDesign location blocks **before** any existing `location /` catch-all:

```bash
sudo nano /etc/nginx/conf.d/shiny.conf    # or vim, etc.
```

Add these two blocks inside the `server { listen 443 ssl ...; }` block, before the existing `location / { ... }`:

```nginx
# ── ProbeDesign Streamlit ─────────────────────────────────
location = /probedesign {
    return 301 /probedesign/;
}

location /probedesign/ {
    proxy_pass         http://probedesign;
    proxy_http_version 1.1;
    proxy_set_header   Upgrade $http_upgrade;
    proxy_set_header   Connection "upgrade";
    proxy_set_header   Host $host;
    proxy_set_header   X-Real-IP $remote_addr;
    proxy_set_header   X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header   X-Forwarded-Proto $scheme;
    proxy_read_timeout 86400;
}
```

Why **before** `location /`: Nginx uses the most specific prefix match. `location /probedesign/` is more specific than `location /`, so it always wins regardless of order — but placing it first makes it visually obvious and avoids any ambiguity.

The `location = /probedesign` (exact match, no slash) issues a redirect to the trailing-slash version. Without this, `location /probedesign/` (prefix match) does **not** match `/probedesign` and nginx falls through to the catch-all, returning a 404 from the Shiny server.

**The complete `shiny.conf` with ProbeDesign added:**

```nginx
server {
    listen 80 default_server;
    server_name shinytest.mvls.gla.ac.uk;
    return 301 https://$host$request_uri;
}

server {
    listen 443 ssl http2 default_server;
    server_name shinytest.mvls.gla.ac.uk;

    ssl_certificate /etc/letsencrypt/live/shinytest.mvls.gla.ac.uk/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/shinytest.mvls.gla.ac.uk/privkey.pem;

    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_prefer_server_ciphers on;
    ssl_ciphers 'TLS_AES_128_GCM_SHA256:TLS_AES_256_GCM_SHA384:TLS_CHACHA20_POLY1305_SHA256:ECDHE+AESGCM:ECDHE+CHACHA20';

    # ── ProbeDesign Streamlit ─────────────────────────────────
    location = /probedesign {
        return 301 /probedesign/;
    }

    location /probedesign/ {
        proxy_pass         http://probedesign;
        proxy_http_version 1.1;
        proxy_set_header   Upgrade $http_upgrade;
        proxy_set_header   Connection "upgrade";
        proxy_set_header   Host $host;
        proxy_set_header   X-Real-IP $remote_addr;
        proxy_set_header   X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header   X-Forwarded-Proto $scheme;
        proxy_read_timeout 86400;
    }

    # ── R Shiny (catch-all) ───────────────────────────────────
    location / {
        proxy_pass http://127.0.0.1:3838;
        proxy_http_version 1.1;
        proxy_redirect off;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        rewrite ^/(.*)$ /$1 break;
    }
}
```

**8e. Allow SELinux to proxy to Streamlit (required on Oracle Linux 9):**

SELinux prevents `httpd` (Nginx) from making outbound connections to non-standard ports by default. Enable permanently:

```bash
sudo setsebool -P httpd_can_network_connect on

# Verify
getsebool httpd_can_network_connect
# Should print: httpd_can_network_connect --> on
```

**8f. Test the config syntax and restart:**

```bash
sudo nginx -t
sudo systemctl restart nginx
```

Expected output from `nginx -t`: `syntax is ok` / `test is successful`. A warning about conflicting `server_name` is acceptable if IT's config and yours both use `_` as a catch-all; what matters is no `[emerg]` errors.

---

### Step 9 — Configure the Firewall

Oracle Linux 9 uses `firewalld`. Since IT handles HTTPS, port 443 may already be open; verify before adding rules:

```bash
sudo firewall-cmd --list-services
# If output already includes: https   → nothing to do
# If not, open it:
sudo firewall-cmd --permanent --add-service=https
sudo firewall-cmd --reload
sudo firewall-cmd --list-services   # should now include: https
```

Port 8501 stays localhost-only (Nginx reaches it internally — never expose it externally).

**University network security group**: On this VM, IT have already configured port 443 access at the network level. If deploying on a different VM, contact IT to confirm HTTPS is allowed inbound.

---

### Step 10 — Verify the Deployment

From another machine on the university network, open a browser and navigate to:

```
https://shinytest.mvls.gla.ac.uk/probedesign/
```

- HTTP (`http://...`) automatically redirects to HTTPS (IT's config)
- `/probedesign` without a trailing slash redirects to `/probedesign/` (our `location = /probedesign` block)
- The R Shiny server continues to work on all other paths (IT's `location /` catch-all is unaffected)

**Test a probe design run:**

1. Paste or upload a FASTA sequence (e.g., from `tests/regression/CDKN1A_32/CDKN1A.fa`)
2. Set probes to 32
3. Enable Pseudogene mask and Genome mask
4. Click **Design Probes**
5. Verify results appear and files are downloadable

**Monitor memory during genome masking** (in a separate SSH session):

```bash
watch -n 1 free -m
```

During genome masking, expect RAM + swap to be heavily used (~3.5 GB for human genome bowtie queries). The existing 3.9 GB swap should handle this, but observe to confirm the process does not get OOM-killed.

---

## Maintenance & Updates

### Updating the App (git pull)

When new features are pushed to the repository:

```bash
# SSH into the VM
cd ~/Github/smFISHProbeDesign

# Pull latest changes
git pull origin master

# If environment.yml changed (new dependencies), update the env:
micromamba activate probedesign
micromamba env update -n probedesign -f environment.yml --prune -y

# If pyproject.toml changed, reinstall the package:
pip install -e .

# Restart the service to pick up changes
sudo systemctl restart probedesign

# Verify
sudo systemctl status probedesign
```

**Quick update** (if only Python code changed, no new dependencies):

```bash
cd ~/Github/smFISHProbeDesign
git pull origin master
sudo systemctl restart probedesign
```

### Viewing Logs

```bash
# Live log stream
sudo journalctl -u probedesign -f

# Last 100 lines
sudo journalctl -u probedesign --no-pager -n 100

# Logs from today only
sudo journalctl -u probedesign --since today

# Search for errors
sudo journalctl -u probedesign --no-pager | grep -i error
```

### Restarting the Service

```bash
sudo systemctl restart probedesign      # Restart
sudo systemctl stop probedesign         # Stop
sudo systemctl start probedesign        # Start
sudo systemctl status probedesign       # Check status
```

### Monitoring Memory

```bash
# One-shot check
free -m

# Continuous monitoring (update every 2 seconds)
watch -n 2 free -m

# Check if OOM killer has acted
sudo dmesg | grep -i "out of memory"
sudo journalctl --since today | grep -i "oom"
```

---

## Troubleshooting

### Service fails to start

```bash
# Check the logs
sudo journalctl -u probedesign --no-pager -n 50

# Common issues:
# 1. Wrong path to streamlit → fix ExecStart in service file
# 2. Wrong username → fix User/Group and paths in service file
# 3. Missing bowtie → check PATH environment in service file
```

After editing the service file:

```bash
sudo systemctl daemon-reload
sudo systemctl restart probedesign
```

### "bowtie: command not found" in service

The systemd service runs without shell initialization, so conda's PATH isn't set automatically. The `Environment="PATH=..."` line in the service file must include the conda env's bin directory. Verify:

```bash
# Find where bowtie is installed
micromamba activate probedesign
which bowtie
# Example: /home/<username>/micromamba/envs/probedesign/bin/bowtie

# The PATH in the service file must include that directory:
# Environment="PATH=/home/<username>/micromamba/envs/probedesign/bin:..."
```

### App accessible on VM but not from other machines

Work through these in order:

**1. Test from inside the VM first:**
```bash
curl -s -o /dev/null -w "%{http_code}\n" http://localhost/probedesign/
```
- `301` or `200` → Nginx routing works locally; the issue is external connectivity or HTTPS
- `502` → Nginx routing works but can't reach Streamlit (check SELinux, check service is running)
- `404` → Nginx not routing to our location block (check config)

**2. Check if the response is an HTTP→HTTPS redirect:**
```bash
curl -v http://localhost/probedesign/ 2>&1 | grep "Location:"
```
If `Location: https://...` appears, IT's config is redirecting all HTTP to HTTPS. Everything must be configured on the HTTPS (port 443) server block — a port-80-only config will never be reached.

**3. Check for IT's existing nginx config:**
```bash
ls /etc/nginx/conf.d/
sudo grep -rn "return 301\|listen 443\|ssl_certificate" /etc/nginx/conf.d/
```
If IT has a file with `return 301 https://` on port 80 and a `listen 443 ssl` block, your location blocks must go inside that HTTPS server block (see Step 8d).

**4. Check for conflicting `server_name` or `default_server`:**
```bash
sudo nginx -T 2>/dev/null | grep -n "server_name\|listen\|default_server"
```
`[warn] conflicting server name` is harmless. `[emerg] duplicate default server` means two server blocks both have `default_server` on the same port — remove one.

**5. SELinux blocking Nginx → Streamlit connection:**
```bash
getsebool httpd_can_network_connect
# If "off": sudo setsebool -P httpd_can_network_connect on
```

**6. Nginx status and logs:**
```bash
sudo systemctl status nginx
sudo tail -20 /var/log/nginx/error.log
```

### Genome masking causes OOM kill

Symptoms: The app hangs or crashes during genome masking. Check:

```bash
sudo dmesg | tail -20   # Look for "Out of memory: Killed process"
```

If this happens, the 3.4 GB RAM + 3.9 GB swap may not be enough under heavy load. Options:

1. **Add more swap** (easiest):
   ```bash
   sudo fallocate -l 4G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile

   # Make persistent across reboots
   echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab

   # Verify
   free -h   # Should now show ~7.9 GB swap total
   ```

2. **Limit concurrent users**: Only one genome masking job should run at a time. Streamlit's default behavior (one Python process per browser tab) makes this unlikely unless multiple users run genome masking simultaneously.

3. **Request more RAM from IT** if genome masking is a critical feature.

### Streamlit shows "Please wait..." forever

```bash
# Check if the process is running
sudo systemctl status probedesign

# Check if port is listening
ss -tlnp | grep 8501

# Try accessing locally on the VM
curl -s http://localhost:8501 | head -5
```

### Updating fails (git pull conflicts)

If someone edited files directly on the VM:

```bash
cd ~/Github/smFISHProbeDesign

# Option 1: Stash local changes, pull, then reapply
git stash
git pull origin master
git stash pop

# Option 2: Discard local changes entirely (destructive)
git checkout -- .
git pull origin master
```

---

## Optional: Install RepeatMasker (Later)

If you later decide to enable automatic repeat masking, follow these steps. This requires ~65 GB free disk space during installation (~56 GB permanent).

### 1. Install RepeatMasker via conda

```bash
micromamba activate probedesign
micromamba install -c bioconda -c conda-forge repeatmasker -y
```

### 2. Download the Dfam database (Mammalia partition)

This is the large download — ~8.9 GB compressed, ~56 GB extracted.

```bash
# Find the RepeatMasker Libraries directory
RMPATH=$(dirname $(which RepeatMasker))/../share/RepeatMasker/Libraries/famdb
echo "RepeatMasker famdb dir: $RMPATH"

cd "$RMPATH"

# Download Partition 0 (root, always needed) — if not already present
# It's usually included with the conda install, but verify:
ls dfam39_full.0.h5 2>/dev/null || {
    curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.0.h5.gz
    gunzip dfam39_full.0.h5.gz
}

# Download Partition 7 (Mammalia — needed for human, mouse, rat)
curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.7.h5.gz
gunzip dfam39_full.7.h5.gz    # This creates a ~56 GB file
```

### 3. Verify RepeatMasker

```bash
micromamba activate probedesign
RepeatMasker -v       # Should print version
RepeatMasker -species human -dir /tmp/rm_test tests/regression/CDKN1A_32/CDKN1A.fa
ls /tmp/rm_test/      # Should contain CDKN1A.fa.masked
```

### 4. Restart the service

```bash
sudo systemctl restart probedesign
```

RepeatMasker will now be detected automatically by the app's prerequisite check. The "Auto RepeatMasker" option will become available in the sidebar.

---

## Quick Reference

| Task | Command |
|------|---------|
| SSH in | `ssh <username>@<vm-ip>` |
| Start screen session | `screen -S probedesign_setup` |
| Detach from screen | `Ctrl-A` then `D` |
| Reattach to screen | `screen -r probedesign_setup` |
| Check service | `sudo systemctl status probedesign` |
| Restart service | `sudo systemctl restart probedesign` |
| View live logs | `sudo journalctl -u probedesign -f` |
| Update app | `cd ~/Github/smFISHProbeDesign && git pull origin master && sudo systemctl restart probedesign` |
| Update env | `micromamba activate probedesign && micromamba env update -n probedesign -f environment.yml --prune -y && pip install -e . && sudo systemctl restart probedesign` |
| Check memory | `free -m` |
| Check Nginx status | `sudo systemctl status nginx` |
| Reload Nginx config | `sudo nginx -t && sudo systemctl restart nginx` |
| Check firewall | `sudo firewall-cmd --list-services` |
| Verify index dir | `ls ~/Github/smFISHProbeDesign/bowtie_indexes/*.1.bt2 \| wc -l` |
| App URL | `https://shinytest.mvls.gla.ac.uk/probedesign/` |
