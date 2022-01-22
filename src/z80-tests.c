
#include "../z80.c"

// Testing code.

#define MEM_SIZE 0x10000

static int load_file(const char* filename, u16 addr) {
  FILE* f = fopen(filename, "r");
  fseek(f, 0, SEEK_END);
  size_t s = ftell(f);
  rewind(f);

  if (s + addr >= MEM_SIZE) {
    fprintf(stderr, "%s too big to fit in memory.\n", filename);
    R 1;
  }

  fread(m + addr, sizeof(u8), s, f);
  fclose(f);
  R 0;
}

u8 port_in(u8 port) { R 0; }

void port_out(u8 port, u8 val) {
  if(val == 1) {
    // Service 2: print character.
    if (c == 2) {
      printf("%c", e);
    }

    // Service 9: print a '$'-terminated string.
    if (c == 9) {
      u16 addr = (d << 8) | e;
      do
        putchar(r8(addr++));
      while (r8(addr) != '$');
    }
  }
}

static void run_test(const char* filename) {
  init();
  memset(m, 0, MEM_SIZE);

  if (load_file(filename, 0x100) != 0)
    R;
  puts(filename);

  pc = 0x100;

  // write 1 to port 0 - handle BDOS services.
  // ld a, 1
  // out 0, a
  // ret
  m[5] = 0x3E;
  m[6] = 0x01;
  m[7] = 0xD3;
  m[8] = 0x00;
  m[9] = 0xC9;

  do step(); while(pc);
  putchar('\n');
}

int main(void) {
  m = malloc(MEM_SIZE);

  run_test("roms/zexdoc.com");
  run_test("roms/zexall.com");
}
