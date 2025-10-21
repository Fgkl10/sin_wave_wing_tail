import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev
from shapely.geometry import LineString

# ---------- Parameters and function definitions (same as original) ----------
M = 0.04
P = 0.3
T = 0.19
a0 = 0.2969
a1 = -0.126
a2 = -0.3516
a3 = 0.2843
a4 = -0.1036
d = 0.3
n = 17

def y_c(x):
    x = np.array(x)
    y = np.zeros_like(x)
    mask1 = (0 <= x) & (x <= P)
    mask2 = (P < x) & (x <= 1)
    y[mask1] = (M/(P**2))*(2*P*x[mask1] - x[mask1]**2)
    y[mask2] = (M/((1-P)**2))*(((1-2*P) + 2*P*x[mask2] - x[mask2]**2))
    return y

def y_dc(x):
    x = np.array(x)
    y = np.zeros_like(x)
    mask1 = (0 <= x) & (x <= P)
    mask2 = (P < x) & (x <= 1)
    y[mask1] = (2*M/(P**2))*(P - x[mask1])
    y[mask2] = (2*M/((1-P)**2))*(P - x[mask2])
    return y

def theta(x):
    return np.arctan(y_dc(x))

def y_t(x):
    return (T/0.2)*(a0*np.sqrt(x) + a1*x + a2*x**2 + a3*x**3 + a4*x**4)

def x_u(x):
    return x - y_t(x)*np.sin(theta(x))

def y_u(x):
    return y_c(x) + y_t(x)*np.cos(theta(x))

def y_l(x):
    return y_c(x) - y_t(x)*np.cos(theta(x))

def A(x):
    return (y_u(x) - y_l(x))/2

b = np.pi / d * n

scale = 150.0
scale_x = scale  # x-axis scaling
scale_y = 130    # y-axis scaling

num_points = 150
x = np.linspace(d, 1, num_points)
y = y_c(x) + A(x) * np.sin(b * x)

x_scaled = x * scale_x
y_scaled = y * scale_y

# ---------- Thicken along the normal direction, generate closed curve ----------
dx = np.gradient(x_scaled)
dy = np.gradient(y_scaled)
length = np.sqrt(dx**2 + dy**2)
tx = dx / length
ty = dy / length
nx = -ty
ny = tx

width = 1

x_outer = x_scaled + width/2 * nx
y_outer = y_scaled + width/2 * ny
x_inner = x_scaled - width/2 * nx
y_inner = y_scaled - width/2 * ny

x_closed = np.concatenate([x_outer, x_inner[::-1]])
y_closed = np.concatenate([y_outer, y_inner[::-1]])


# ---------- Self-intersection check function ----------
def is_self_intersect(x, y):
    line = LineString(np.column_stack([x, y]))
    return not line.is_simple

# ---------- Automatically enlarge s, find lower and upper bounds, then binary search for minimal non-self-intersecting s ----------
def find_min_s(x_closed, y_closed, s_start=1.0, step_factor=10, max_s=1e5, tol=1e-2):
    s_lower = s_start
    s_upper = s_start
    # 1. Incrementally search for non-self-intersecting interval
    while s_upper <= max_s:
        tck, u = splprep([x_closed, y_closed], s=s_upper, per=True)
        unew = np.linspace(0, 1, len(x_closed))
        x_bspline, y_bspline = splev(unew, tck)
        if not is_self_intersect(x_bspline, y_bspline):
            break
        s_lower = s_upper
        s_upper *= step_factor
    else:
        raise RuntimeError("No non-self-intersecting s found; please increase max_s or adjust curve parameters.")

    print(f"Incremental initial search interval: s_lower={s_lower}, s_upper={s_upper}")

    # 2. Binary search for minimal non-self-intersecting s
    while s_upper - s_lower > tol:
        s_mid = (s_lower + s_upper) / 2
        tck, u = splprep([x_closed, y_closed], s=s_mid, per=True)
        unew = np.linspace(0, 1, len(x_closed))
        x_bspline, y_bspline = splev(unew, tck)
        if is_self_intersect(x_bspline, y_bspline):
            s_lower = s_mid
        else:
            s_upper = s_mid

    # s_upper is the minimal non-self-intersecting s
    tck, u = splprep([x_closed, y_closed], s=s_upper, per=True)
    unew = np.linspace(0, 1, len(x_closed))
    x_bspline, y_bspline = splev(unew, tck)
    return s_upper, x_bspline, y_bspline

# ---------- Main procedure ----------
s_start = 1      # Initial s parameter
step_factor = 10   # Step multiplier
max_s = 100000000    # Maximum s parameter
tol = 0.1        # Binary search convergence tolerance

min_s, x_bspline, y_bspline = find_min_s(x_closed, y_closed, s_start, step_factor, max_s, tol)
print(f"Minimal non-self-intersecting s parameter: {min_s:.4f}")

# ---------- Add 0.1 automatically when writing ----------
min_s_out = min_s + 0.1

# ---------- Closed curve check and automatic closure ----------
# Check if the distance between the first and last point is less than threshold (e.g., 0.1); if not, automatically close the curve
close_threshold = 0.1
dist_closure = np.hypot(x_bspline[0] - x_bspline[-1], y_bspline[0] - y_bspline[-1])
if dist_closure > close_threshold:
    # Automatically add the first point as the last point to close the curve
    x_bspline = np.append(x_bspline, x_bspline[0])
    y_bspline = np.append(y_bspline, y_bspline[0])
    print(f"Curve is not closed, closed automatically.")
else:
    print(f"Curve is already closed, no action needed.")

# ---------- Output xyz file (n_width_numpoints_s format, skip points closer than 0.1, keep two decimals) ----------
out_fname = f"{n}_{width}_{num_points}_{min_s_out:.1f}.txt"
written = 0
skipped = 0
last_x = None
last_y = None
with open(out_fname, "w") as f:
    for xi, yi in zip(x_bspline, y_bspline):
        # Distance check
        if last_x is not None and last_y is not None:
            dist = np.hypot(xi - last_x, yi - last_y)
            if dist < 0.1:
                skipped += 1
                continue
        # Write point
        f.write(f"{xi:.2f},{yi:.2f},0.00\n")
        written += 1
        last_x, last_y = xi, yi

    # Add one more line to close the curve (repeat the starting point)
    f.write(f"{x_bspline[0]:.2f},{y_bspline[0]:.2f},0.00\n")
    written += 1

print(f"File written: {out_fname}, lines written: {written} (including closure point), skipped close points: {skipped}")

# ---------- Visualization ----------
plt.figure(figsize=(10, 6))
plt.fill(x_bspline, y_bspline, color='skyblue', alpha=0.7, label="B-spline Smoothed Curve")
plt.plot(x_scaled, y_scaled, color='r', lw=1, label="Original Center Line")
plt.xlabel("x")
plt.ylabel("y")
plt.title(f"Fixed Thickness Wavy Curve (B-spline Smoothed, Min s={min_s:.2f})")
plt.legend()
plt.axis("equal")
plt.grid(True)
plt.show()