
#include "../z80.c"
#include <termios.h>
#include <sys/ioctl.h>

// Testing code.

#define MEM_SIZE 0x10000

void erm() {
  struct termios term;
  tcgetattr(0, &term);
  term.c_lflag &= ~(ICANON | ECHO);
  tcsetattr(0, TCSANOW, &term);
}

void drm() {
  struct termios term;
  tcgetattr(0, &term);
  term.c_lflag |= ICANON | ECHO;
  tcsetattr(0, TCSANOW, &term);
}

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

// Stubs.
u8 port_in(u8 port) { R 0; }
void port_out(u8 port, u8 val) { }

int basic_started = 0;

// Initialisation routine.
void rst00() {
  sp = 0x20ED;
  printf("BasicBIOS for Microsoft BASIC v4.7 by Kamila Szewczyk.\n");
  if(!basic_started)
    goto cold_start;
  printf("Cold/warm boot? (C/W) ");
  int c = getchar();
  if(c == 'W')
    goto warm_start;

  cold_start:
    basic_started = 1;
    pc = 0x0150;
    return;
  warm_start:
    pc = 0x0153;
    return;
}

int crlf = 0;

_Bool kbhit() {
    int bw;
    ioctl(0, FIONREAD, &bw);
    return bw > 0;
}

void rst08() { putchar(a); ret(); }
void rst10() {
  if(crlf) {
    crlf = 0;
    a = '\n';
    ret(); return;
  }
  a = getchar();
  if(a == '\n')
    a = '\r', crlf = 1;
  if(a == 127) { a = '\b'; printf("\b "); }
  ret();
}
void rst18() { zf = !kbhit(); ret(); }
void rst38() { ret(); } // interrupt stub.

static void run_basic(const char* filename) {
  init();
  memset(m, 0, MEM_SIZE);

  if (load_file(filename, 0x150) != 0)
    R;
  puts(filename);

  pc = 0x0000;

  do {
    switch(pc){
    case 0x0000: rst00(); break; // rst 0x00
    case 0x0008: rst08(); break; // rst 0x08
    case 0x0010: rst10(); break; // rst 0x10
    case 0x0018: rst18(); break; // rst 0x18
    case 0x0038: rst38(); break; // rst 0x38
    default: step();
    }
  } while(pc);
  putchar('\n');
}

int main(void) {
  m = malloc(MEM_SIZE);
  atexit(drm);
  erm();

  run_basic("roms/msbasic.com");
}
