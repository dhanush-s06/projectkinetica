"""
demo.py — Demonstrations of the orbital.py library.

Generates three interactive HTML visualizations:
  1. earth_moon.html  — Earth–Moon system (one lunar month)
  2. solar.html       — Inner solar system (Sun + 4 planets, two Earth-years)
  3. binary.html      — A binary star system with a third orbiting body
"""

from orbital import Body, System, inner_solar_system

# ---------------------------------------------------------------------------
# 1. Earth–Moon system
# ---------------------------------------------------------------------------
def earth_moon():
    s = System(name="Earth & Moon")
    earth = Body("Earth", mass=5.972e24, position=[0, 0], velocity=[0, 0],
                 color="#5b8def", radius=14)
    moon = Body("Moon", mass=7.347e22,
                position=[3.844e8, 0], velocity=[0, 1022.0],
                color="#d8d3c4", radius=5)
    s.add(earth, moon)
    # Counter the Moon's pull on Earth so the system stays roughly centered
    earth.velocity[1] = -moon.velocity[1] * moon.mass / earth.mass
    s.simulate(duration=30 * 86_400, dt=600, sample_every=4)  # ~30 days
    s.export_html("earth_moon.html", title="Earth & Moon")
    print("✓ earth_moon.html")

# ---------------------------------------------------------------------------
# 2. Inner solar system
# ---------------------------------------------------------------------------
def solar():
    s = inner_solar_system()
    # 2 Earth years, 12-hour steps, sample every 2 steps -> 1-day resolution
    s.simulate(duration=2 * 365.25 * 86_400, dt=12 * 3600, sample_every=2)
    s.export_html("solar.html", title="The Inner Solar System")
    print("✓ solar.html")

# ---------------------------------------------------------------------------
# 3. Binary star + planet (a chaotic three-body dance)
# ---------------------------------------------------------------------------
def binary():
    s = System(name="Binary Star with Companion")
    M = 2.0e30
    sep = 5e10
    # Two equal stars in mutual circular orbit around their barycenter
    v_star = (6.6743e-11 * M / (4 * sep)) ** 0.5
    s.add(Body("Castor", M, [-sep, 0], [0,  v_star], color="#e8b87c", radius=14))
    s.add(Body("Pollux", M, [ sep, 0], [0, -v_star], color="#cd6a4a", radius=14))
    # Companion in a much wider orbit; treat the binary as a 2M point mass.
    r_comp = 4e11
    v_comp = (6.6743e-11 * (2 * M) / r_comp) ** 0.5  # circular velocity
    s.add(Body("Companion", 1e25, [0, r_comp], [-v_comp, 0],
               color="#5b8def", radius=6))
    s.simulate(duration=6 * 365.25 * 86_400, dt=3600 * 6, sample_every=2)
    s.export_html("binary.html", title="Binary Star & Companion")
    print("✓ binary.html")


if __name__ == "__main__":
    earth_moon()
    solar()
    binary()
    print("\nDone. Open any .html in a browser.")
