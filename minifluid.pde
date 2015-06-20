FluidGrid fluidGrid;

void setup() {
  size(1280, 720, P2D);
  
  fluidGrid = new FluidGrid(10, 1280, 720);
}
void draw() {
  // physics
  fluidGrid.solve();
  
  background(255);
  fluidGrid.draw();
}
