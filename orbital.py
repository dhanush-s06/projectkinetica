"""
orbital.py — a small, dependency-light orbital mechanics library that exports
self-contained, browser-native interactive visualizations.

Usage
-----
    from orbital import Body, System

    sys = System(name="Earth–Moon")
    earth = Body("Earth", mass=5.972e24, position=[0, 0], velocity=[0, 0],
                 color="#5b8def", radius=14)
    moon  = Body("Moon",  mass=7.347e22, position=[3.844e8, 0], velocity=[0, 1022],
                 color="#d8d3c4", radius=5)
    sys.add(earth, moon)
    sys.simulate(duration=30 * 86_400, dt=600, sample_every=4)
    sys.export_html("earth_moon.html")

The simulation uses a 4th-order Runge–Kutta integrator and supports any
number of bodies in 2D. Output is a single HTML file with no external
dependencies (other than two web fonts).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence

# Gravitational constant (m^3 kg^-1 s^-2)
G = 6.67430e-11


# ---------------------------------------------------------------------------
# Bodies & Systems
# ---------------------------------------------------------------------------

@dataclass
class Body:
    """A point-mass gravitating body in 2D."""
    name: str
    mass: float
    position: Sequence[float]
    velocity: Sequence[float]
    color: str = "#e8b87c"
    radius: float = 8.0  # display radius in pixels (visual only)
    trail: List[List[float]] = field(default_factory=list)

    def __post_init__(self):
        self.position = [float(self.position[0]), float(self.position[1])]
        self.velocity = [float(self.velocity[0]), float(self.velocity[1])]


class System:
    """A gravitational N-body system with an RK4 integrator."""

    def __init__(self, name: str = "Orbital System"):
        self.name = name
        self.bodies: List[Body] = []
        self.time: float = 0.0
        self._sample_dt: float = 0.0  # seconds between recorded samples

    # ----- construction ----------------------------------------------------

    def add(self, *bodies: Body) -> "System":
        for b in bodies:
            self.bodies.append(b)
        return self

    # ----- physics ---------------------------------------------------------

    def _accelerations(self, positions: List[List[float]]) -> List[List[float]]:
        """Compute gravitational accelerations on each body (pure-Python; small N)."""
        n = len(positions)
        acc = [[0.0, 0.0] for _ in range(n)]
        for i in range(n):
            xi, yi = positions[i]
            for j in range(n):
                if i == j:
                    continue
                dx = positions[j][0] - xi
                dy = positions[j][1] - yi
                r2 = dx * dx + dy * dy
                if r2 == 0.0:
                    continue
                inv_r3 = 1.0 / (r2 * math.sqrt(r2))
                m_j = self.bodies[j].mass
                acc[i][0] += G * m_j * dx * inv_r3
                acc[i][1] += G * m_j * dy * inv_r3
        return acc

    def _rk4_step(self, dt: float) -> None:
        n = len(self.bodies)
        x = [list(b.position) for b in self.bodies]
        v = [list(b.velocity) for b in self.bodies]

        def add(a, b, s):
            return [[a[i][0] + s * b[i][0], a[i][1] + s * b[i][1]] for i in range(n)]

        a1 = self._accelerations(x)
        k1x = [list(vi) for vi in v]
        k1v = a1

        a2 = self._accelerations(add(x, k1x, 0.5 * dt))
        k2x = add(v, k1v, 0.5 * dt)
        k2v = a2

        a3 = self._accelerations(add(x, k2x, 0.5 * dt))
        k3x = add(v, k2v, 0.5 * dt)
        k3v = a3

        a4 = self._accelerations(add(x, k3x, dt))
        k4x = add(v, k3v, dt)
        k4v = a4

        for i, body in enumerate(self.bodies):
            body.position[0] += dt / 6.0 * (k1x[i][0] + 2 * k2x[i][0] + 2 * k3x[i][0] + k4x[i][0])
            body.position[1] += dt / 6.0 * (k1x[i][1] + 2 * k2x[i][1] + 2 * k3x[i][1] + k4x[i][1])
            body.velocity[0] += dt / 6.0 * (k1v[i][0] + 2 * k2v[i][0] + 2 * k3v[i][0] + k4v[i][0])
            body.velocity[1] += dt / 6.0 * (k1v[i][1] + 2 * k2v[i][1] + 2 * k3v[i][1] + k4v[i][1])

        self.time += dt

    # ----- simulation ------------------------------------------------------

    def simulate(self, duration: float, dt: float, sample_every: int = 1) -> "System":
        """Integrate forward by `duration` seconds in steps of `dt`,
        recording each body's position every `sample_every` steps."""
        if dt <= 0:
            raise ValueError("dt must be positive")
        if sample_every < 1:
            raise ValueError("sample_every must be >= 1")
        steps = int(round(duration / dt))
        # Record initial state
        for b in self.bodies:
            b.trail.append([b.position[0], b.position[1]])
        for i in range(1, steps + 1):
            self._rk4_step(dt)
            if i % sample_every == 0:
                for b in self.bodies:
                    b.trail.append([b.position[0], b.position[1]])
        self._sample_dt = dt * sample_every
        return self

    # ----- conveniences ----------------------------------------------------

    @staticmethod
    def circular_orbit_velocity(central_mass: float, radius: float) -> float:
        """Speed required for a circular orbit around a fixed mass at given radius."""
        return math.sqrt(G * central_mass / radius)

    # ----- export ----------------------------------------------------------

    def export_html(self, path: str, *, title: str | None = None,
                    show_trails: bool = True, show_grid: bool = True) -> str:
        """Write a self-contained HTML file with the simulation embedded.
        Returns the absolute path to the written file."""
        title = title or self.name
        # Build the data payload
        bodies_payload = [{
            "name": b.name,
            "mass": b.mass,
            "color": b.color,
            "radius": b.radius,
            "trail": b.trail,
        } for b in self.bodies]
        payload = {
            "title": title,
            "sample_dt": self._sample_dt,
            "bodies": bodies_payload,
            "options": {
                "show_trails": bool(show_trails),
                "show_grid": bool(show_grid),
            },
        }
        # Load and inject into the template
        template_path = Path(__file__).with_name("template.html")
        if not template_path.exists():
            raise FileNotFoundError(
                f"template.html not found next to orbital.py at {template_path}"
            )
        html = template_path.read_text(encoding="utf-8")
        html = html.replace("/*__PAYLOAD__*/null", json.dumps(payload))
        html = html.replace("__TITLE__", title)
        out = Path(path).resolve()
        out.write_text(html, encoding="utf-8")
        return str(out)


# ---------------------------------------------------------------------------
# Convenience factory: build well-known systems quickly
# ---------------------------------------------------------------------------

def inner_solar_system() -> System:
    """A pre-baked inner solar system (Sun + 4 inner planets), 2D approximation."""
    AU = 1.495978707e11
    sun_mass = 1.98892e30
    s = System(name="Inner Solar System")
    s.add(Body("Sun", sun_mass, [0, 0], [0, 0], color="#ffce6e", radius=18))
    for name, dist_au, mass, color, r in [
        ("Mercury", 0.387, 3.301e23, "#b8a89a", 4),
        ("Venus",   0.723, 4.867e24, "#e8c992", 6),
        ("Earth",   1.000, 5.972e24, "#5b8def", 7),
        ("Mars",    1.524, 6.417e23, "#cd6a4a", 5),
    ]:
        r_m = dist_au * AU
        v = System.circular_orbit_velocity(sun_mass, r_m)
        s.add(Body(name, mass, [r_m, 0], [0, v], color=color, radius=r))
    return s
