#!/bin/bash
# Re-download every paper listed in papers/README.md into this directory.
# Existing files are kept. Run from anywhere: bash papers/fetch.sh
set -euo pipefail
cd "$(dirname "$0")"
PAPERS=(
  0905.2280 1604.07818 1610.03816 1702.08016 1811.02584 1904.01911
  1911.10171 1912.00004 2006.13983 2012.12326 2203.00700 2206.01083
  2302.11882 2303.00792 2310.16994 2312.15598 2401.08555 2403.01950 2403.05980
  2404.03535 2406.11078 2407.09396 2407.12988 2409.11324 2501.09423
  2505.08116 2506.18054v3 2512.07715 2601.09801 2607.05678v1
  2607.13665v1 2607.14215v1 2607.14215v2
)
for id in "${PAPERS[@]}"; do
  [ -f "$id.pdf" ] && continue
  echo "fetching $id"
  curl -sfL -o "$id.pdf" "https://arxiv.org/pdf/$id"
  sleep 2   # arXiv rate limit
done
# The target paper's LaTeX source (line-number citations in the code use it)
if [ ! -d arXiv-2607.05678v1 ]; then
  curl -sfL https://arxiv.org/e-print/2607.05678v1 -o src.tar.gz
  mkdir -p arXiv-2607.05678v1 && tar -xzf src.tar.gz -C arXiv-2607.05678v1
  rm src.tar.gz
fi
