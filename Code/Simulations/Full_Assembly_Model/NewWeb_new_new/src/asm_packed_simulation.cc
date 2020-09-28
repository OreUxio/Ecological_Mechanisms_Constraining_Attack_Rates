#include "packed_simulation.h"

"Sorry, this won't work. We need to convert everything into an .s file,
and just assemble it using 'as'.  Let's get this function below bound
as a simple c-function, to avoid naming problems.  We need the class
membership nonly for postinc, which is easily done in assembler"

double packed_simulation::
asm_fast_sandwich_product(const double * v1_in,const double * v2_in,double factor,
			  char * __restrict__ tsk){
  double last;
  double vec1_in[bs];
  register double * vec1 asm("r9");
  char * matrix_end = postinc<char *>(tsk);  // strangely, we need this!!
  register aligned_value_t * atsk asm ("rcx") = aligned_value_pointer(tsk);
  register double sum asm ("xmm0") = 0;
  register const double * v1 asm("rsi") = v2_in;
  register const double * v2 asm("rdx") = v1_in;
  do{
  asm volatile
    ( ".intel_syntax noprefix\n"
      "lea    r11,[r9+0x10]\n"
      "lea    r10,[r9+0x20]\n"
"2:	xor eax,eax\n"
      "jmp    1\n"
      "nop\n"
      "nop\n"
      "nop\n"
      "nop\n"
"1:	movsx  r8,WORD PTR [rcx+rax*1+0x2]\n"
      "movsx  rdi,WORD PTR [rcx+rax*1]\n"
      "movsd  xmm1,QWORD PTR [rdx+r8*8]\n"
      "mulsd  xmm1,QWORD PTR [rsi+rdi*8]\n"
      "movsd  QWORD PTR [r9+rax*2],xmm1\n"
      "add    rax,0x4\n"
      "cmp    rax,0x20\n"
      "jne    1\n"
      "movapd xmm2,XMMWORD PTR [rcx+0x20]\n"
      "movapd xmm1,XMMWORD PTR [rcx+0x30]\n"
      "mulpd  xmm2,XMMWORD PTR [r9]\n"
      "mulpd  xmm1,XMMWORD PTR [r11]\n"
      "addpd  xmm2,xmm1\n"
      "movapd xmm1,XMMWORD PTR [rcx+0x40]\n"
      "mulpd  xmm1,XMMWORD PTR [r10]\n"
      "addpd  xmm2,xmm1\n"
      "movapd xmm1,XMMWORD PTR [rcx+0x50]\n"
      "add    rcx,0x60\n"
      "mulpd  xmm1,XMMWORD PTR [r10+0x10]\n"
      "addpd  xmm1,xmm2\n"
      "haddpd xmm1,xmm1\n"
      "addsd  xmm0,xmm1\n" : : : 
      "xmm1", "xmm2", "xmm3", "eax", "r8", "r10", "r11", "rdi");
      // "movsd  xmm1,QWORD PTR [rcx-0x8]\n"
      // "mulsd  xmm1,xmm3\n"
      // "comisd xmm1,xmm0\n"
      // "ja     .L2"
    last = *(((value_t *)atsk)-1);
  }while(last*factor>sum);
  return sum;
}
