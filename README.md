# Ornament of Pressure

Real-time Islamic geometric ornament visualizer with audio reactivity.

A native Windows C++ application that renders intricate Islamic geometric patterns — reacting to live microphone input. The ornaments are based on constructions from Daud Sutton's *Islamic Design: A Genius for Geometry* (Wooden Books, 2007), a concise guide to the geometric principles underlying centuries of Islamic art and architecture.

![Screenshot](screenshot.png)

## About the Patterns

Islamic geometric art is built from a small set of principles — regular polygon tessellations, star polygons, and interlace — that generate astonishing complexity. This project implements six patterns drawn from Sutton's book:

| Preset | Source | Construction |
|--------|--------|-------------|
| **Star & Hexagon** | pp. 2–5 | Seal of Solomon — overlapping triangles on a hexagonal grid |
| **Breath of the Compassionate** | pp. 8–9 | Expanding/contracting double squares on a square grid |
| **Eight-Fold Rosettes** | pp. 10–11 | Octagonal rosettes with kite-shaped petals |
| **Six of One (12-fold)** | pp. 16–19 | 12-pointed stars from dodecagon edge midpoints |
| **Umm al-Girih (10-fold)** | pp. 34–37 | Decagonal stars with pentagonal sub-patterns |
| **Zillij (8-fold)** | pp. 26–29 | Octagram seal {8/3} with nested octagons |

Each pattern is generated procedurally using Hankin/PIC construction on regular grids, then traced into closed contours for rendering and filling.

## Features

- **6 book-authentic ornament presets** with labeled sidebar thumbnails
- **Audio-reactive geometry** — FFT analysis of microphone input drives pattern depth, contact angle, and rotation speed (bass/mids/highs)
- **Geometric morphing** — switching ornaments interpolates vertex positions between patterns (no fade, pure geometry)
- **Filled regions** — half of all closed contours are filled with semi-transparent white using stencil even-odd rule
- **Adjustable sensitivity** — scroll wheel or arrow keys control audio reactivity gain
- **Closed contours only** — all geometry consists of proper closed loops (no dangling lines)

## Controls

| Input | Action |
|-------|--------|
| Click thumbnail / Keys 1-6 | Switch ornament |
| Scroll wheel / Up-Down arrows | Adjust audio sensitivity |
| Escape | Quit |

## Building

Requires MSYS2 with ucrt64 toolchain:

```bash
pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-SDL2 mingw-w64-ucrt-x86_64-glew

export PATH="/c/msys64/ucrt64/bin:$PATH"
g++ -O2 -std=c++17 main.cpp -o ornament.exe -lmingw32 -lSDL2main -lSDL2 -lglew32 -lopengl32
```

Make sure `SDL2.dll` and `glew32.dll` are in the same directory as the executable. A prebuilt `ornament.exe` is included for quick demos.

## Architecture

Single-file C++ (`main.cpp`, ~1300 lines):

- **Geometry** — Hankin construction on square, triangular, and hexagonal grids; star polygons {n/k}; polar star boundaries
- **Audio** — SDL2 audio capture, 1024-sample Cooley-Tukey FFT with Hann windowing, 3-band analysis (bass/mids/highs)
- **Rendering** — OpenGL 2.1 fixed-function pipeline, VBO + glMultiDrawArrays, stencil-based polygon fills
- **Chain building** — spatial hash graph, degree-1 pruning, straightest-path closed loop tracing
- **Morphing** — chain resampling and vertex-position interpolation with smoothstep easing
- **Text** — minimal stroke font for sidebar labels (no font library dependency)

## References

- Daud Sutton, *Islamic Design: A Genius for Geometry*, Wooden Books, 2007 — [Internet Archive](https://archive.org/details/islamicdesigngen0000sutt)

## Web Version

`index.html` contains an earlier web-based version with basic Hankin construction on a triangle grid. Open directly in a browser.
