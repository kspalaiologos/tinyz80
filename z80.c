
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

// ---------------------------------------
// ABBREVIATIONS

typedef int8_t i8;
typedef uint8_t u8;
typedef uint16_t u16;
typedef _Bool bv;
#define SI static inline
#define R return

// Defining a single-expression function body.
// Equivalent declarations:
// u8 f(u8 a)_(a * 2)
// u8 f(u8 a) { return a * 2; }
#define _(a...) {return({a;});}

// ---------------------------------------
// PROCESSOR STATE.

// The z80 memory.
u8 * m = NULL;

// Instruction counter, stack pointer, index registers.
u16 pc, sp, ix, iy;

// "WZ" register.
// https://www.grimware.org/lib/exe/fetch.php/documentations/devices/z80/z80.memptr.eng.txt
u16 wz;

// Main and exchange registers.
u8 a, b, c, d, e, h, l;
u8 exa, exb, exc, exd, exe, exh, exl, exf;

// Interrupt vectors and the refresh register.
// Note that the refresh register is useless but since
// it can be loaded to a, we must implement it properly.
// Some programs use it as a semi-random entropy source.
u8 i, r;

// Flags: sign, zero, halfcarry, parity, negative, carry.
// yf and xf are undocumented flags. In this emulator,
// they are hardwired to 3rd and 5th bits of the result.
// http://www.z80.info/z80sflag.htm
bv sf, zf, yf, hf, xf, pf, nf, cf;

// Interrupt flip-flops. `iff_set' is set when `EI' is executed.
// https://eduinf.waw.pl/inf/retro/004_z80_inst/0005.php
u8 iff_set;
bv iff1, iff2;

// The interrupt mode. Set by `IM', is either 0, 1 or 2.
// Mode 2 is rarely used, but it needs to be implemented.
u8 im;
u8 int_vec;
bv int_pending, nmi_pending;

bv halted;

// ---------------------------------------
// API

// Port I/O functions supplied by the user.
extern u8 port_in(u8 port);
extern void port_out(u8 port, u8 val);

// Exposed API.
void init();
void step();
void gen_nmi();
void gen_int(u8 data);

// ---------------------------------------
// MEMORY / PROCESSOR STATE OPERATIONS.

// Extracting n-th bit of a number, byte-wise memory access.
#define bit(n, val) (((val) >> (n)) & 1)
#define r8(a) m[a]
#define w8(a, v) m[a] = v

// Higher/lower byte of a u16.
#define HI(a) ((a) >> 8)
#define LO(a) ((a) & 0xFF)

// Memory access with word granularity.
SI u16 r16(u16 a) _ ((r8(a + 1) << 8) | r8(a))
SI void w16(u16 a, u16 v) { w8(a, LO(v)); w8(a + 1, HI(v)); }

// Stack operations.
SI void psh(u16 v) { w16(sp -= 2, v); }
SI u16 pop() _ (r16((sp += 2) - 2))

// Querying program data and incrementing the instruction pointer.
SI u8 p8() _ (r8(pc++))
SI u16 p16() _ (r16((pc += 2) - 2))

SI void set_bc(u16 v) { b = HI(v); c = LO(v); }
SI void set_hl(u16 v) { h = HI(v); l = LO(v); }
SI void set_de(u16 v) { d = HI(v); e = LO(v); }

// Setting / getting the "full" flags register.
// Order: cf nf pf xf hf yf zf sf.
SI u8 get_f() _ (
  cf | (nf << 1) | (pf << 2) | (xf << 3)
| (hf << 4) | (yf << 5) | (zf << 6) | (sf << 7))

SI void set_f(u8 val) {
  cf = bit(0, val);
  nf = bit(1, val);
  pf = bit(2, val);
  xf = bit(3, val);
  hf = bit(4, val);
  yf = bit(5, val);
  zf = bit(6, val);
  sf = bit(7, val);
}

// 16-bit register pairs.
SI u16 get_bc() _ ((b << 8) | c)
SI u16 get_hl() _ ((h << 8) | l)
SI u16 get_de() _ ((d << 8) | e)
SI u16 get_af() _ ((a << 8) | get_f())

// ---------------------------------------
// INSTRUCTION IMPLEMENTATIONS.

// increments R, preserving the highest bit.
SI void inc_r() { r = (r & 0x80) | ((r + 1) & 0x7f); }

// find the carry bit in addition on bits n and n-1
SI bv carry(int n, u16 a, u16 b, bv cy) _ (bit(n, (a + b + cy) ^ a ^ b))

// Check whether the parity of a byte is even.
SI bv parity(u8 x) {
  x ^= x >> 4;
  x ^= x >> 2;
  x ^= x >> 1;
  return (~x) & 1;
}

// Branching
SI void jmp(u16 a) { wz = pc = a; }
SI void call(u16 dest) { psh(pc); jmp(dest); }
SI void cjmp(bv c) { u16 t = p16(); if (c) jmp(t); wz = t; }
SI void ccall(bv c) { u16 t = p16(); if (c) call(t); wz = t; }
SI void ret() { jmp(pop()); }
SI void cret(bv c) { if (c) ret(); }

// Relative jumps.
SI void jr(i8 disp) { wz = pc += disp; }
SI void cjr(bv c) { if (c) jr(p8()); else pc++; }

// Updating XF and YF in two ways.
SI void xyf1(u8 x) { xf = bit(3, x); yf = bit(5, x); }
SI void xyf2(u8 x) { xf = bit(3, x); yf = bit(1, x); }

// Adjust S and Z flags.
SI void szf8(u8 x) { zf = x == 0; sf = x >> 7; }
SI void szf16(u16 x) { zf = x == 0; sf = x >> 15; }

// Arithmetics
SI u8 add8(u8 a, u8 b, bv cy) {
  u8 res = a + b + cy;
  hf = carry(4, a, b, cy);
  pf = carry(7, a, b, cy) != carry(8, a, b, cy);
  cf = carry(8, a, b, cy);
  nf = 0;
  szf8(res);
  xyf1(res);
  R res;
}

// Subtraction by addition of negated value.
SI u8 sub8(u8 a, u8 b, bv cy) {
  u8 v = add8(a, ~b, !cy);
  cf = !cf; hf = !hf;
  nf = 1;
  R v;
}

SI u16 add16(u16 a, u16 b, bv cy) {
  u8 lo = add8(a, b, cy);
  
  u16 res = (add8(HI(a), HI(b), cf) << 8) | lo;
  zf = res == 0; wz = a + 1; R res;
}

SI u16 sub16(u16 a, u16 b, bv cy) {
  u8 lo = sub8(a, b, cy);
  
  u16 res = (sub8(HI(a), HI(b), cf) << 8) | lo;
  zf = res == 0; wz = a + 1; R res;
}

// Add a word to hl.
SI void addhl(u16 val) {
  bv nsf = sf; bv nzf = zf; bv npf = pf;
  set_hl(add16(get_hl(), val, 0));
  sf = nsf; zf = nzf; pf = npf;
}

// Add a word to a native 16-bit register.
SI void addiz(u16 * reg, u16 val) {
  bv nsf = sf; bv nzf = zf; bv npf = pf;
  *reg = add16(*reg, val, 0);
  sf = nsf; zf = nzf; pf = npf;
}

// Add/subtract to HL with carry.
SI void adchl(u16 v) { u16 q = add16(get_hl(), v, cf); szf16(q); set_hl(q); }
SI void sbchl(u16 v) { u16 q = sub16(get_hl(), v, cf); szf16(q); set_hl(q); }
SI u8 inc(u8 a) { bv ncf = cf; u8 q = add8(a, 1, 0); cf = ncf; R q; }
SI u8 dec(u8 a) { bv ncf = cf; u8 q = sub8(a, 1, 0); cf = ncf; R q; }

// Bit operations.
#define bit_op(name, op, halfcarry) \
  SI void name(u8 val) { \
    u8 res = a op val; \
    szf8(res); \
    xyf1(res); \
    hf = halfcarry; \
    pf = parity(res); \
    nf = cf = 0; \
    a = res; \
  }

bit_op(land, &, 1)
bit_op(lor,  |, 0)
bit_op(lxor, ^, 0)

#undef bit_op

// compare with A.
SI void cmpa(u8 val) { sub8(a, val, 0); xyf1(val); }

// A stub for 0xCB opcs.
#define op_cbh(name, op...) \
  SI u8 name(u8 val) { \
    op \
    szf8(val); \
    xyf1(val); \
    pf = parity(val); \
    nf = hf = 0; \
    R val; \
  }

// Rotate left/right with carry.
op_cbh(rlc, bv old = val >> 7; val = (val << 1) | old; cf = old;)
op_cbh(rrc, bv old = val & 1; val = (val >> 1) | (old << 7); cf = old;)

// Simple shifts.
op_cbh(rl, bv ncf = cf; cf = val >> 7; val = (val << 1) | ncf;)
op_cbh(rr, bv nc = cf; cf = val & 1; val = (val >> 1) | (nc << 7);)

// Shift and preserve sign (possibly: set the bit to 1).
op_cbh(sla, cf = val >> 7; val <<= 1;)
op_cbh(sll, cf = val >> 7; val <<= 1; val |= 1;)
op_cbh(sra, cf = val & 1; val = (val >> 1) | (val & 0x80);)
op_cbh(srl, cf = val & 1; val >>= 1;)

// Test n-th bit.
SI u8 bt(u8 val, u8 n) {
  u8 res = val & (1 << n);
  szf8(res); xyf1(val);
  hf = 1; pf = zf; nf = 0; R res;
}

// Altering 16-bit register pairs.
#define dbc(n) set_bc(get_bc() + (n))
#define dhl(n) set_hl(get_hl() + (n))
#define dde(n) set_de(get_de() + (n))

// Copying data.
SI void ldi() {
  u8 val = r8(get_hl());
  
  w8(get_de(), val);
  
  dhl(1); dde(1); dbc(-1);
  
  xyf2(val + a);
  
  nf = hf = 0;
  pf = get_bc() > 0;
}

SI void ldd() { ldi(); dhl(-2); dde(-2); }

// Comparison
SI void cpi() {
  bv ncf = cf;
  u8 v = sub8(a, r8(get_hl()), 0);
  dhl(1); dbc(-1);
  xyf2(v - hf);
  pf = get_bc() != 0;
  cf = ncf;
  wz += 1;
}

SI void cpd() { cpi(); dhl(-2); wz -= 2; }

// PIO
SI void inr(u8 * r) {
  *r = port_in(c);
  szf8(*r);
  pf = parity(*r);
  nf = hf = 0;
}

SI void adji() {
  dhl(1);
  zf = --b == 0;
  nf = 1;
  wz = get_bc() + 1;
}

SI void ini() { w8(get_hl(), port_in(c)); adji(); }
SI void outi() { port_out(c, r8(get_hl())); adji(); }

SI void ind() { ini(); dhl(-2); wz = get_bc() - 2; }
SI void outd() { outi(); dhl(-2); wz = get_bc() - 2; }

// Decimal adjust.
SI void daa() {
  u8 bcd = 0;
  
  if ((a & 0x0F) > 0x09 || hf)
    bcd += 0x06;
  
  if (a > 0x99 || cf)
    bcd += 0x60, cf = 1;
  
  if (nf) {
    hf = hf && (a & 0x0F) < 0x06;
    bcd = -bcd;
  } else {
    hf = (a & 0x0F) > 0x09;
  }
  
  pf = parity(a += bcd);
  xyf1(a); szf8(a);
}

// Displacement computation. Updates the WZ pair.
SI u16 dp(u16 b, i8 d) _ (wz = b + d)

// ---------------------------------------
// EXECUTION.

SI void exec(u8 opc);
SI void exec_cb(u8 opc);
SI void exec_cb2(u8 opc, u16 addr);
SI void exec_ed(u8 opc);
SI void exec_ind(u8 opc, u16 * ir);

// User interface.
void init() {
  pc = ix = iy = wz = 0;
  sp = 0xFFFF;
  
  // AF = 0xFFFF, zero the rest.
  a = 0xFF;
  sf = zf = 1;
  xf = yf = 1;
  hf = pf = 1;
  nf = cf = 1;
  
  b = c = d = e = h = l = 0;
  exa = exb = exc = exd = exe = exh = exl = 0;
  
  i = r = 0;
  
  iff_set = 0;
  im = 0;
  iff1 = iff2 = 0;
  halted = 0;
  int_pending = nmi_pending = 0;
  int_vec = 0;
}

// Execute an instruction and handle interrupts.
void step() {
  exec(halted ? 0x00 : p8());
  
  if (iff_set) {
    iff_set = 0;
    iff1 = iff2 = 1;
    R;
  }
  
  if (nmi_pending) {
    nmi_pending = halted = iff1 = 0;
    inc_r(); call(0x66); R;
  }
  
  if (int_pending && iff1) {
    int_pending = halted = 0;
    iff1 = iff2 = 0;
    inc_r();
    
    // Dispatch the interrupt based on the interrupt mode.
    switch (im) {
    case 0: exec(int_vec); R;
    case 1: call(0x38); R;
    case 2: call(r16((i << 8) | int_vec)); R;
    }
  }
}

// Schedule a NMI
void gen_nmi() {
  nmi_pending = 1;
}

// Schedule an interrupt.
void gen_int(u8 data) {
  int_pending = 1;
  int_vec = data;
}

// ---------------------------------------
// INSTRUCTION DECODERS

#define H(n, c...) case n: c; break;

// Adjust WZ after operation.
SI void awz1(u16 val) { wz = (a << 8) | LO(val + 1); }

// Execute a regular opcode.
SI void exec(u8 opc) {
  u16 t1; bv t2; u8 t3;
  
  inc_r();
  
  switch (opc) {
  // ld X, Y block.
  H(0x47, b=a) H(0x40, b=b) H(0x41, b=c) H(0x42, b=d) H(0x43, b=e) H(0x44, b=h) H(0x45, b=l)
  H(0x57, d=a) H(0x50, d=b) H(0x51, d=c) H(0x52, d=d) H(0x53, d=e) H(0x54, d=h) H(0x55, d=l)
  H(0x67, h=a) H(0x60, h=b) H(0x61, h=c) H(0x62, h=d) H(0x63, h=e) H(0x64, h=h) H(0x65, h=l)
  H(0x4F, c=a) H(0x48, c=b) H(0x49, c=c) H(0x4A, c=d) H(0x4B, c=e) H(0x4C, c=h) H(0x4D, c=l)
  H(0x5F, e=a) H(0x58, e=b) H(0x59, e=c) H(0x5A, e=d) H(0x5B, e=e) H(0x5C, e=h) H(0x5D, e=l)
  H(0x6F, l=a) H(0x68, l=b) H(0x69, l=c) H(0x6A, l=d) H(0x6B, l=e) H(0x6C, l=h) H(0x6D, l=l)
  H(0x7F, a=a) H(0x78, a=b) H(0x79, a=c) H(0x7A, a=d) H(0x7B, a=e) H(0x7C, a=h) H(0x7D, a=l)
  
  // ld X, imm
  H(0x0E, c=p8()) H(0x06, b=p8())
  H(0x1E, e=p8()) H(0x16, d=p8())
  H(0x2E, l=p8()) H(0x26, h=p8())
  H(0x3E, a=p8()) H(0x36, w8(get_hl(),p8()))
  
  // ld X, (HL) and ld (HL), X
  H(0x6E, l=r8(get_hl()))
  H(0x7E, a=r8(get_hl())) H(0x46, b=r8(get_hl()))
  H(0x4E, c=r8(get_hl())) H(0x56, d=r8(get_hl()))
  H(0x5E, e=r8(get_hl())) H(0x66, h=r8(get_hl()))
  
  H(0x77, w8(get_hl(),a)) H(0x70, w8(get_hl(),b))
  H(0x71, w8(get_hl(),c)) H(0x72, w8(get_hl(),d))
  H(0x73, w8(get_hl(),e)) H(0x74, w8(get_hl(),h))
  H(0x75, w8(get_hl(),l))
  
  // ld (bc), a / ld (de), a / ld *word, a
  H(0x02, w8(get_bc(),a);awz1(get_bc()))
  H(0x12, w8(get_de(),a);awz1(get_de()))
  H(0x32, t1=p16();w8(t1,a);awz1(t1))
  
  // ld a, (bc) / ld a, (de) / ld a, *word
  H(0x0A, a=r8(get_bc());wz=get_bc()+1)
  H(0x1A, a=r8(get_de());wz=get_de()+1)
  H(0x3A, wz=p16()+1;a=r8(wz-1))
  
  // ld r16, imm word
  H(0x01, set_bc(p16())) H(0x11, set_de(p16()))
  H(0x21, set_hl(p16())) H(0x31, sp=p16())
  
  // ld r16, *word
  H(0x2A, set_hl(r16(t1=p16()));wz=t1+1)
  H(0x22, w16(t1=p16(),get_hl());wz=t1+1)
  
  // ld sp, hl
  H(0xF9, sp=get_hl())
  
  // ex de/(sp), hl
  H(0xEB, t1=get_de();set_de(get_hl());set_hl(t1))
  H(0xE3, t1=r16(sp);w16(sp,get_hl());set_hl(wz=t1))
  
  // add block
  #define op(opc, src) case opc: a = add8(a, src, 0); break;
  op(0x87, a)
  op(0x80, b)
  op(0x81, c)
  op(0x82, d)
  op(0x83, e)
  op(0x84, h)
  op(0x85, l)
  op(0x86, r8(get_hl()))
  op(0xC6, p8())
  #undef op
  
  // add with carry
  #define op(opc, src) case opc: a = add8(a, src, cf); break;
  op(0x8F, a)
  op(0x88, b)
  op(0x89, c)
  op(0x8A, d)
  op(0x8B, e)
  op(0x8C, h)
  op(0x8D, l)
  op(0x8E, r8(get_hl()))
  op(0xCE, p8())
  #undef op
  
  // sub block
  #define op(opc, src) case opc: a = sub8(a, src, 0); break;
  op(0x97, a)
  op(0x90, b)
  op(0x91, c)
  op(0x92, d)
  op(0x93, e)
  op(0x94, h)
  op(0x95, l)
  op(0x96, r8(get_hl()))
  op(0xD6, p8())
  #undef op
  
  // sub with carry
  #define op(opc, src) case opc: a = sub8(a, src, cf); break;
  op(0x9F, a)
  op(0x98, b)
  op(0x99, c)
  op(0x9A, d)
  op(0x9B, e)
  op(0x9C, h)
  op(0x9D, l)
  op(0x9E, r8(get_hl()))
  op(0xDE, p8())
  #undef op
  
  // add hl, ...
  H(0x09, addhl(get_bc())) H(0x19, addhl(get_de()))
  H(0x29, addhl(get_hl())) H(0x39, addhl(sp))
  
  // di, ei, nop, halt
  H(0xF3, iff1=iff2=0)
  H(0xFB, iff_set=1)
  H(0x00, )
  H(0x76, halted=1)
  
  // inc/dec block.
  H(0x3C, a=inc(a)) H(0x04, b=inc(b))
  H(0x0C, c=inc(c)) H(0x14, d=inc(d))
  H(0x1C, e=inc(e)) H(0x24, h=inc(h))
  H(0x2C, l=inc(l))
  H(0x34, w8(get_hl(),inc(r8(get_hl()))))
  
  H(0x3D, a=dec(a)) H(0x05, b=dec(b))
  H(0x0D, c=dec(c)) H(0x15, d=dec(d))
  H(0x1D, e=dec(e)) H(0x25, h=dec(h))
  H(0x2D, l=dec(l))
  H(0x35, w8(get_hl(),dec(r8(get_hl()))))
  
  // inc/dec r16
  H(0x03, dbc(1)) H(0x13, dde(1)) H(0x23, dhl(1)) H(0x33, ++sp)
  H(0x0B, dbc(-1)) H(0x1B, dde(-1)) H(0x2B, dhl(-1)) H(0x3B, --sp)
  
  // decimal adjust
  H(0x27, daa())
  
  // cpl/scf/ccf
  H(0x2F, a=~a;nf=hf=1;xyf1(a))
  H(0x37, cf=1;nf=hf=0;xyf1(a))
  H(0x3F, hf=cf;cf=!cf;nf=0;xyf1(a))
  
  // rlca, rrca, rrl
  H(0x07, cf=a>>7;a=(a<<1)|cf;nf=hf=0;xyf1(a))
  H(0x17, t2=cf;cf=a>>7;a=(a<<1)|t2;nf=hf=0;xyf1(a))
  H(0x0F, cf=a&1;a=(a>>1)|(cf<<7);nf=hf=0;xyf1(a))
  H(0x1F, t2=cf;cf=a&1;a=(a>>1)|(t2<<7);nf=hf=0;xyf1(a))
  
  // bit operations
  #define op(opA, opB, opC, opD, opE, opH, opL, opDHL, opI, name) \
    case opA: name(a); break; \
    case opB: name(b); break; \
    case opC: name(c); break; \
    case opD: name(d); break; \
    case opE: name(e); break; \
    case opH: name(h); break; \
    case opL: name(l); break; \
    case opDHL: name(r8(get_hl())); break; \
    case opI: name(p8()); break;
  op(0xA7, 0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xE6, land)
  op(0xB7, 0xB0, 0xB1, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6, 0xF6, lor)
  op(0xAF, 0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAD, 0xAE, 0xEE, lxor)
  op(0xBF, 0xB8, 0xB9, 0xBA, 0xBB, 0xBC, 0xBD, 0xBE, 0xFE, cmpa)
  #undef op
  
  // jmp imm
  H(0xC3, jmp(p16()))
  
  // conditional jumps.
  H(0xC2, cjmp(!zf)) H(0xCA, cjmp(zf))
  H(0xD2, cjmp(!cf)) H(0xDA, cjmp(cf))
  H(0xE2, cjmp(!pf)) H(0xEA, cjmp(pf))
  H(0xF2, cjmp(!sf)) H(0xFA, cjmp(sf))
  
  H(0x10, cjr(--b)) H(0x18, pc+=(i8)p8())
  H(0x20, cjr(!zf)) H(0x28, cjr(zf))
  H(0x30, cjr(!cf)) H(0x38, cjr(cf))
  
  H(0xE9, pc=get_hl()) H(0xCD, call(p16()))
  
  // call/ret
  H(0xC4, ccall(!zf)) H(0xCC, ccall(zf))
  H(0xD4, ccall(!cf)) H(0xDC, ccall(cf))
  H(0xE4, ccall(!pf)) H(0xEC, ccall(pf))
  H(0xF4, ccall(!sf)) H(0xFC, ccall(sf))
  
  H(0xC9, ret())
  
  H(0xC0, cret(!zf)) H(0xC8, cret(zf))
  H(0xD0, cret(!cf)) H(0xD8, cret(cf))
  H(0xE0, cret(!pf)) H(0xE8, cret(pf))
  H(0xF0, cret(!sf)) H(0xF8, cret(sf))
  
  // rst - 1 byte calls to fixed address.
  // similar to 8086's `int', except we
  // call the IVT, not the address in it.
  H(0xC7, call(0x00)) H(0xCF, call(0x08))
  H(0xD7, call(0x10)) H(0xDF, call(0x18))
  H(0xE7, call(0x20)) H(0xEF, call(0x28))
  H(0xF7, call(0x30)) H(0xFF, call(0x38))
  
  // psh
  H(0xC5, psh(get_bc())) H(0xD5, psh(get_de()))
  H(0xE5, psh(get_hl())) H(0xF5, psh(get_af()))
  
  // pop
  H(0xC1, set_bc(pop()))
  H(0xD1, set_de(pop()))
  H(0xE1, set_hl(pop()))
  H(0xF1, t1=pop();a=t1>>8;set_f(t1&0xFF))
  
  // in/out
  H(0xDB, t3=a;a=port_in(p8());wz=(t3<<8)|(a+1))
  H(0xD3, t3=p8();port_out(t3,a);wz=(t3+1)|(a<<8))
  
  // swapping exhange/main registers.
  case 0x08: {
    u8 na = a, nf = get_f();
    a = exa; set_f(exf);
    exa = na; exf = nf;
    break;
  }
  
  case 0xD9: {
    u8 nb = b, nc = c, nd = d, ne = e, nh = h, nl = l;
    
    b = exb;
    c = exc;
    d = exd;
    e = exe;
    h = exh;
    l = exl;
    
    exb = nb;
    exc = nc;
    exd = nd;
    exe = ne;
    exh = nh;
    exl = nl;
    
    break;
  }
  
  H(0xCB, exec_cb(p8())) H(0xDD, exec_ind(p8(), &ix))
  H(0xED, exec_ed(p8())) H(0xFD, exec_ind(p8(), &iy))
  
  default: fprintf(stderr, "Unknown opcode: 0x%02X\n", opc);
  }
}

// Execute an opc operating on IX or IY.
SI void exec_ind(u8 opc, u16 * ir) {
  u16 t1;
  
  inc_r();
  
  #define IHI HI(*ir)
  #define ILO LO(*ir)
  #define IDP dp(*ir, p8())
  
  switch (opc) {
  // Stack operations.
  H(0xE1, *ir = pop()) H(0xE5, psh(*ir))
  
  // Jumps
  H(0xE9, jmp(*ir))
  
  // Arithmetics.
  H(0x09, addiz(ir, get_bc())) H(0x19, addiz(ir, get_de()))
  H(0x29, addiz(ir, *ir)) H(0x39, addiz(ir, sp))
  
  // hi/lo math.
  // add/adc a, IHI/ILO
  H(0x84, a = add8(a, IHI, 0)) H(0x85, a = add8(a, ILO, 0))
  H(0x8C, a = add8(a, IHI, cf)) H(0x8D, a = add8(a, ILO, cf))
  
  // add/adc/sub/sbc a, byte *(ir + imm)
  H(0x86, a = add8(a, r8(IDP), 0)) H(0x8E, a = add8(a, r8(IDP), cf))
  H(0x96, a = sub8(a, r8(IDP), 0)) H(0x9E, a = sub8(a, r8(IDP), cf))
  
  // sub/sbc a, IHI/ILO
  H(0x94, a = sub8(a, IHI, 0)) H(0x95, a = sub8(a, ILO, 0))
  H(0x9C, a = sub8(a, IHI, cf)) H(0x9D, a = sub8(a, ILO, cf))
  
  #define op3(opA, opB, opC, f) H(opA, f(r8(IDP))) H(opB, f(IHI)) H(opC, f(ILO))
  op3(0xA6, 0xA4, 0xA5, land)
  op3(0xAE, 0xAC, 0xAD, lxor)
  op3(0xB6, 0xB4, 0xB5, lor)
  op3(0xBE, 0xBC, 0xBD, cmpa)
  #undef op3
  
  H(0x23, ++*ir) H(0x2B, --*ir)
  
  // inc/dec *(ix/iy + imm)
  H(0x34, t1=IDP;w8(t1,inc(r8(t1))))
  H(0x35, t1=IDP;w8(t1,dec(r8(t1))))
  
  // inc/dec IHI/ILO
  H(0x24, *ir=ILO|(inc(IHI)<<8))
  H(0x25, *ir=ILO|(dec(IHI)<<8))
  H(0x2C, *ir=(IHI<<8)|inc(ILO))
  H(0x2D, *ir=(IHI<<8)|dec(ILO))
  
  // loading IX/IY.
  H(0x2A, *ir=r16(p16()))
  H(0x22, w16(p16(),*ir))
  H(0x21, *ir=p16())
  
  // Load imm to (ix/iy + imm)
  H(0x36, t1=IDP;w8(t1,p8()))
  
  // ld (ix/iy + imm), r8
  H(0x70, w8(IDP,b)) H(0x71, w8(IDP,c))
  H(0x72, w8(IDP,d)) H(0x73, w8(IDP,e))
  H(0x74, w8(IDP,h)) H(0x75, w8(IDP,l))
  H(0x77, w8(IDP,a))
  
  // loading to registers.
  #define op(ob, oc, od, oe, oh, ol, oa, var) \
    H(ob, b = var) H(oc, c = var) H(od, d = var) H(oe, e = var) \
    H(oh, h = var) H(ol, l = var) H(oa, a = var)
  #define op2(ob, oc, od, oe, oa, var) \
    H(ob, b = var) H(oc, c = var) H(od, d = var) H(oe, e = var) \
    H(oa, a = var)
  op(0x46, 0x4E, 0x56, 0x5E, 0x66, 0x6E, 0x7E, r8(IDP))
  op2(0x44, 0x4C, 0x54, 0x5C, 0x7C, IHI)
  op2(0x45, 0x4D, 0x55, 0x5D, 0x7D, ILO)
  #undef op
  #undef op2
  
  // ld IHI/ILO r/imm/IHI/ILO
  H(0x67, *ir=ILO|(a<<8))
  H(0x60, *ir=ILO|(b<<8))
  H(0x61, *ir=ILO|(c<<8))
  H(0x62, *ir=ILO|(d<<8))
  H(0x63, *ir=ILO|(e<<8))
  H(0x26, *ir=ILO|(p8()<<8))
  H(0x64, )
  H(0x65, *ir=(ILO<<8)|ILO)
  
  H(0x6F, *ir=(IHI<<8)|a)
  H(0x68, *ir=(IHI<<8)|b)
  H(0x69, *ir=(IHI<<8)|c)
  H(0x6A, *ir=(IHI<<8)|d)
  H(0x6B, *ir=(IHI<<8)|e)
  H(0x2E, *ir=(IHI<<8)|p8())
  H(0x6D, )
  H(0x6C, *ir=(IHI<<8)|IHI)
  
  // ld sp, ix/iy
  H(0xF9, sp=*ir)
  
  // ex (sp), ix/iy
  H(0xE3, t1=r16(sp);w16(sp,*ir);wz=*ir=t1)
  
  // two byte prefix
  H(0xCB, t1=IDP;exec_cb2(p8(),t1))
  
  default:
    // non-prefixed opc & decrement r
    exec(opc);
    r = (r & 0x80) | ((r - 1) & 0x7f);
  }
  
  #undef IDP
  #undef ILO
  #undef IHI
}

// executes a CB opcode
SI void exec_cb(u8 opc) {
  inc_r();
  
  // https://eduinf.waw.pl/inf/retro/004_z80_inst/0014.php
  u8 dk = (opc >> 6) & 0b11; // op kind
  u8 da = (opc >> 3) & 0b111; // auxiliary / op0 kind type
  u8 dr = opc & 0b111; // data
  
  // auxiliary storage for data under hl.
  u8 hlr = 0;
  
  // pointer to the storage.
  u8* reg = 0;
  switch (dr) {
  H(0,reg=&b) H(1,reg=&c) H(2,reg=&d) H(3,reg=&e) H(4,reg=&h)
  H(5,reg=&l) H(6, hlr=r8(get_hl());reg=&hlr) H(7,reg=&a)
  }
  
  switch (dk) {
  // Rotation
  case 0:
    switch (da) {
    #define op(n,f) case n: *reg = f(*reg); break;
    op(0, rlc) op(1, rrc)
    op(2, rl) op(3, rr)
    op(4, sla) op(5, sra)
    op(6, sll) op(7, srl)
    #undef op
    }
    break;
  // Testing bits
  case 1:
    bt(*reg, da);
    
    // Odd edge case. See the WZ comment.
    if (dr == 6)
      xyf1(HI(wz));
    break;
  H(2, *reg&=~(1<<da))
  H(3, *reg|=1<<da)
  }
  
  // Write back to hl ptr if needed.
  if (reg == &hlr)
    w8(get_hl(), hlr);
}

// Execute a CB opcode with regards to ix/iy.
SI void exec_cb2(u8 opc, u16 addr) {
  u8 val = r8(addr);
  u8 res = 0;
  
  u8 dk = (opc >> 6) & 0b11; // op kind
  u8 da = (opc >> 3) & 0b111; // auxiliary / op0 kind type
  u8 dr = opc & 0b111; // data
  
  switch (dk) {
  // Rotation
  case 0:
    switch (da) {
    #define op(n,f) case n: res = f(val); break;
    op(0, rlc) op(1, rrc)
    op(2, rl) op(3, rr)
    op(4, sla) op(5, sra)
    op(6, sll) op(7, srl)
    #undef op
    }
    break;
  case 1:
    res = bt(val, da);
    xyf1(HI(addr));
    break;
  H(2, res=val&~(1<<da))
  H(3, res=val|(1<<da))
  
  default: fprintf(stderr, "Invalid IX/IY CB-prefixed opcode: %02X\n", opc);
  }
  
  if (dk != 1 && dr != 6)
    switch (dr) {
    H(0, b=res) H(1, c=res) H(2, d=res) H(3, e=res)
    H(4, h=res) H(5, l=res) H(6, w8(get_hl(), res))
    H(7, a=res)
    }
  
  if (dk != 1)
    w8(addr, res);
}

SI void rot_ep() {
  nf = hf = 0;
  xyf1(a);
  szf8(a);
  pf = parity(a);
  wz = get_hl() + 1;
}

// Execute a ED opcode.
SI void exec_ed(u8 opc) {
  u8 t1; u16 t2;
  inc_r();
  switch (opc) {
  H(0x47, i=a)
  H(0x4F, r=a)
  
  // ld a, i/r
  H(0x57, a=i;szf8(a);hf=nf=0;pf=iff2)
  H(0x5F, a=r;szf8(a);hf=nf=0;pf=iff2)
  
  // retn
  case 0x45: case 0x55: case 0x65: case 0x75:
  case 0x5D: case 0x6D: case 0x7D:
    iff1 = iff2; ret(); break;
  // reti
  H(0x4D, ret())
  
  // ldi & ldir, ldd & lddr
  H(0xA0, ldi())
  H(0xB0, ldi();if(get_bc())wz=--pc,--pc;)
  H(0xA8, ldd())
  H(0xB8, ldd();if(get_bc())wz=--pc,--pc;)
  
  // cpi & cpd
  H(0xA1, cpi()) H(0xA9, cpd())
  
  // cpir & cpdr
  H(0xB1, cpi();if(get_bc()&&!zf)wz=--pc,--pc;else++wz)
  H(0xB9, cpd();if(get_bc()&&!zf)pc-=2;else++wz)
  
  // in r, (c)
  H(0x40, inr(&b)) H(0x48, inr(&c))
  H(0x50, inr(&d)) H(0x58, inr(&e))
  H(0x60, inr(&h)) H(0x68, inr(&l))
  
  // in + discard
  H(0x70, inr(&t1));
  
  // in a, (c). note different wz behavior
  H(0x78, inr(&a);wz=get_bc()+1)
  
  // ini/ind
  H(0xA2, ini()) H(0xAA, ind())
  
  // inir / indr
  H(0xB2, ini();if(b)pc-=2)
  H(0xBA, ind();if(b)pc-=2)
  
  // out
  H(0x41, port_out(c, b)) H(0x49, port_out(c, c))
  H(0x51, port_out(c, d)) H(0x59, port_out(c, e))
  H(0x61, port_out(c, h)) H(0x69, port_out(c, l))
  H(0x71, port_out(c, 0))
  
  // again, out a is special-cased for wz.
  H(0x79, port_out(c,a);wz=get_bc()+1)
  
  // outi/outd + repeated versions
  H(0xA3, outi()) H(0xAB, outd())
  H(0xB3, outi();if(b)pc-=2)
  H(0xBB, outd();if(b)pc-=2)
  
  // sbc/adc hl, r16
  H(0x42, sbchl(get_bc())) H(0x52, sbchl(get_de()))
  H(0x62, sbchl(get_hl())) H(0x72, sbchl(sp))
  H(0x4A, adchl(get_bc())) H(0x5A, adchl(get_de()))
  H(0x6A, adchl(get_hl())) H(0x7A, adchl(sp))
  
  // ld [imm], r16
  H(0x43, w16(t2=p16(),get_bc());wz=t2+1)
  H(0x53, w16(t2=p16(),get_de());wz=t2+1)
  H(0x63, w16(t2=p16(),get_hl());wz=t2+1)
  H(0x73, w16(t2=p16(),sp);wz=t2+1)
  
  // ld r16, [imm]
  H(0x4B, set_bc(r16(t2=p16()));wz=t2+1)
  H(0x5B, set_de(r16(t2=p16()));wz=t2+1)
  H(0x6B, set_hl(r16(t2=p16()));wz=t2+1)
  H(0x7B, sp = r16(t2=p16());wz=t2+1)
  
  // neg
  case 0x44: case 0x54: case 0x64: case 0x74:
  case 0x4C: case 0x5C: case 0x6C: case 0x7C:
    a = sub8(0, a, 0);
    break;
  
  // im x
  case 0x46: case 0x66: im = 0; break;
  case 0x56: case 0x76: im = 1; break;
  case 0x5E: case 0x7E: im = 2; break;
  
  // rrd / rld
  case 0x67: {
    u8 na = a, val = r8(get_hl());
    a = (na & 0xF0) | (val & 0xF);
    w8(get_hl(), (val >> 4) | (na << 4));
    rot_ep();
    break;
  }
  
  case 0x6F: {
    u8 na = a, val = r8(get_hl());
    a = (na & 0xF0) | (val >> 4);
    w8(get_hl(), (val << 4) | (na & 0xF));
    rot_ep();
    break;
  }
  
  default: fprintf(stderr, "unknown ED opcode: %02X\n", opc); break;
  }
}

#undef bit
#undef H
#undef HI
#undef LO
