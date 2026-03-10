"""Plot per-interval BEAST timing from screen log to diagnose throttling."""

import re
import sys
import matplotlib.pyplot as plt

log_file = sys.argv[1] if len(sys.argv) > 1 else "out/h3n2/beast-log.txt"

iterations = []
minutes = []

with open(log_file) as f:
    for line in f:
        m = re.search(r"(\d+)\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+(\d+)m(\d+)s/Msamples", line)
        if m:
            iterations.append(int(m.group(1)))
            minutes.append(int(m.group(2)) + int(m.group(3)) / 60)

if not iterations:
    print("No timing data found in", log_file)
    sys.exit(1)

# Compute per-interval wall time from cumulative average
# The reported time is cumulative average, so per-interval = n*avg(n) - (n-1)*avg(n-1)
per_interval = []
for i in range(len(minutes)):
    if i == 0:
        per_interval.append(minutes[0])
    else:
        dt = (i + 1) * minutes[i] - i * minutes[i - 1]
        per_interval.append(dt)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

ax1.plot(iterations, minutes, linewidth=0.5)
ax1.set_ylabel("Cumulative avg (min/M)")
ax1.set_title(f"BEAST timing: {log_file}")
ax1.axhline(y=minutes[0], color="green", linestyle="--", alpha=0.5, label=f"Initial: {minutes[0]:.1f} min/M")
ax1.legend()

ax2.plot(iterations, per_interval, linewidth=0.5, color="tab:orange")
ax2.set_ylabel("Per-interval (min/M)")
ax2.set_xlabel("Iteration")
ax2.axhline(y=minutes[0], color="green", linestyle="--", alpha=0.5)

plt.tight_layout()
out_path = log_file.replace(".txt", "-timing.png")
plt.savefig(out_path, dpi=150)
print(f"Saved to {out_path}")
print(f"Initial: {minutes[0]:.1f} min/M, Final: {minutes[-1]:.1f} min/M")
print(f"Per-interval range: {min(per_interval):.1f} - {max(per_interval):.1f} min/M")
