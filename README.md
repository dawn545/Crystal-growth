# CrystalGrowth

This project simulates dendritic crystal growth using the Kobayashi model. It was initially based on the [jklae/CrystalGrowth](https://github.com/jklae/CrystalGrowth) project, which used DXViewer for rendering. I modified it to use FreeGLUT for cross-platform compatibility and improved accessibility.

I hope this simple project will serve as a helpful tutorial for beginners interested in graphics, specifically those exploring ice crystal growth simulations.
If you have any suggestions or questions, feel free to share!

crystal.exe is 3d version
## Features

- **Dendritic Crystal Growth Simulation**: Models crystal growth in a 2D grid using the Kobayashi model, incorporating anisotropic material properties and thermal dynamics.
- **FreeGLUT Integration**: Replaces DXViewer with FreeGLUT for rendering and interaction, ensuring compatibility with a wider range of platforms.
- **Physical Modeling**: Includes calculations for gradient, Laplacian, and anisotropy effects on crystal growth, temperature field evolution, and phase transitions.

## Reference
Kim, Y., & Lin, S. (2003). Visual Simulation of Ice Crystal Growth. Eurographics/SIGGRAPH Symposium on Computer Animation.

## Runtime effect
<img width="1074" height="1105" alt="屏幕截图 2026-02-28 175208" src="https://github.com/user-attachments/assets/d2c17ca3-a653-4713-b939-00327da0b5a4" />
<img width="561" height="496" alt="屏幕截图 2026-03-15 015233" src="https://github.com/user-attachments/assets/b6d753ce-d4d9-4848-8c40-2ee0e3b96edd" />
<img width="601" height="555" alt="屏幕截图 2026-02-26 181034" src="https://github.com/user-attachments/assets/7c77b66a-5548-46de-aaaa-6f30b5becfc4" />

## Compilation Command

To compile the project, use the following command:

```bash
g++ main.cpp Kobayashi.cpp -I. -lopengl32 -lfreeglut -o main.exe
```bash
