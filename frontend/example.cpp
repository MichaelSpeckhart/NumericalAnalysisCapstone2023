/* em++ example.cpp -s WASM=1 -o example.wasm */

extern "C" {
  int add(int a, int b) {
    return a + b;
  }
}
