	.file	"packed_simulation_model.c"
	.section	.debug_abbrev,"",@progbits
.Ldebug_abbrev0:
	.section	.debug_info,"",@progbits
.Ldebug_info0:
	.section	.debug_line,"",@progbits
.Ldebug_line0:
	.text
.Ltext0:
	.p2align 4,,15
.globl aligned_pointer
	.type	aligned_pointer, @function
aligned_pointer:
.LFB12:
	.file 1 "packed_simulation_model.c"
	.loc 1 15 0
	.cfi_startproc
.LVL0:
	.loc 1 15 0
	leaq	15(%rdi), %rax
	andq	$-16, %rax
	.loc 1 19 0
	ret
	.cfi_endproc
.LFE12:
	.size	aligned_pointer, .-aligned_pointer
	.p2align 4,,15
.globl aligned_value_pointer
	.type	aligned_value_pointer, @function
aligned_value_pointer:
.LFB13:
	.loc 1 22 0
	.cfi_startproc
.LVL1:
	.loc 1 22 0
	leaq	15(%rdi), %rax
	andq	$-16, %rax
	.loc 1 25 0
	ret
	.cfi_endproc
.LFE13:
	.size	aligned_value_pointer, .-aligned_value_pointer
	.p2align 4,,15
.globl asm_fast_sandwich_product
	.type	asm_fast_sandwich_product, @function
asm_fast_sandwich_product:
.LFB14:
	.loc 1 32 0
	.cfi_startproc
.LVL2:
	pushq	%rbp
.LCFI0:
	.cfi_def_cfa_offset 16
	movapd	%xmm0, %xmm3
	movq	%rsp, %rbp
	.cfi_offset 6, -16
.LCFI1:
	.cfi_def_cfa_register 6
.LBB5:
.LBB7:
	.loc 1 24 0
	addq	$23, %rdx
.LVL3:
.LBE7:
.LBE5:
	.loc 1 32 0
	subq	$96, %rsp
	.cfi_escape 0x2e,0x50
.LBB9:
.LBB6:
	.loc 1 24 0
	andq	$-16, %rdx
.LBE6:
.LBE9:
	.loc 1 34 0
	leaq	15(%rsp), %r9
	.loc 1 32 0
	movq	%fs:40, %rax
	movq	%rax, -8(%rbp)
	xorl	%eax, %eax
	.loc 1 34 0
	andq	$-16, %r9
.LBB10:
.LBB8:
	.loc 1 24 0
	xorpd	%xmm0, %xmm0
.LVL4:
	leaq	16(%r9), %r11
	leaq	32(%r9), %r10
	.p2align 4,,10
	.p2align 3
.L7:
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L6:
.LBE8:
.LBE10:
.LBB11:
	.loc 1 45 0
	movswq	2(%rdx,%rax),%r8
	movswq	(%rdx,%rax),%rcx
	movsd	(%rsi,%r8,8), %xmm1
	mulsd	(%rdi,%rcx,8), %xmm1
	movsd	%xmm1, (%r9,%rax,2)
	addq	$4, %rax
	.loc 1 41 0
	cmpq	$32, %rax
	jne	.L6
	.loc 1 51 0
	movapd	32(%rdx), %xmm2
	movsd	(%r9), %xmm1
	movsd	(%r11), %xmm4
	movhpd	8(%r9), %xmm1
	movhpd	8(%r11), %xmm4
	mulpd	%xmm1, %xmm2
	movapd	48(%rdx), %xmm1
	mulpd	%xmm4, %xmm1
	movsd	(%r10), %xmm4
	addpd	%xmm1, %xmm2
	movhpd	8(%r10), %xmm4
	movapd	64(%rdx), %xmm1
	mulpd	%xmm4, %xmm1
	movsd	16(%r10), %xmm4
	addpd	%xmm1, %xmm2
	movhpd	24(%r10), %xmm4
	movapd	80(%rdx), %xmm1
	.loc 1 32 0
	addq	$96, %rdx
.LVL5:
	.loc 1 51 0
	mulpd	%xmm4, %xmm1
	addpd	%xmm2, %xmm1
	haddpd	%xmm1, %xmm1
	addsd	%xmm1, %xmm0
.LBE11:
	.loc 1 55 0
	movsd	-8(%rdx), %xmm1
.LVL6:
	mulsd	%xmm3, %xmm1
	comisd	%xmm0, %xmm1
	ja	.L7
	.loc 1 57 0
	movq	-8(%rbp), %rax
	xorq	%fs:40, %rax
	jne	.L13
.LVL7:
	leave
	ret
.LVL8:
.L13:
	.p2align 4,,9
	.p2align 3
	call	__stack_chk_fail
.LVL9:
	.cfi_endproc
.LFE14:
	.size	asm_fast_sandwich_product, .-asm_fast_sandwich_product
.globl bs
	.section	.rodata
	.align 4
	.type	bs, @object
	.size	bs, 4
bs:
	.long	8
	.text
.Letext0:
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
.LLST2:
	.quad	.LFB14-.Ltext0
	.quad	.LCFI0-.Ltext0
	.value	0x2
	.byte	0x77
	.sleb128 8
	.quad	.LCFI0-.Ltext0
	.quad	.LCFI1-.Ltext0
	.value	0x2
	.byte	0x77
	.sleb128 16
	.quad	.LCFI1-.Ltext0
	.quad	.LFE14-.Ltext0
	.value	0x2
	.byte	0x76
	.sleb128 16
	.quad	0x0
	.quad	0x0
.LLST3:
	.quad	.LVL2-.Ltext0
	.quad	.LVL9-.Ltext0
	.value	0x1
	.byte	0x55
	.quad	0x0
	.quad	0x0
.LLST4:
	.quad	.LVL2-.Ltext0
	.quad	.LVL9-.Ltext0
	.value	0x1
	.byte	0x54
	.quad	0x0
	.quad	0x0
.LLST5:
	.quad	.LVL2-.Ltext0
	.quad	.LVL4-.Ltext0
	.value	0x1
	.byte	0x61
	.quad	.LVL4-.Ltext0
	.quad	.LVL9-.Ltext0
	.value	0x1
	.byte	0x64
	.quad	0x0
	.quad	0x0
.LLST6:
	.quad	.LVL2-.Ltext0
	.quad	.LVL3-.Ltext0
	.value	0x1
	.byte	0x51
	.quad	0x0
	.quad	0x0
.LLST7:
	.quad	.LVL4-.Ltext0
	.quad	.LVL7-.Ltext0
	.value	0x1
	.byte	0x61
	.quad	.LVL8-.Ltext0
	.quad	.LVL9-.Ltext0
	.value	0x1
	.byte	0x61
	.quad	0x0
	.quad	0x0
.LLST8:
	.quad	.LVL3-.Ltext0
	.quad	.LVL5-.Ltext0
	.value	0x1
	.byte	0x51
	.quad	.LVL6-.Ltext0
	.quad	.LVL9-.Ltext0
	.value	0x1
	.byte	0x51
	.quad	0x0
	.quad	0x0
	.section	.debug_info
	.long	0x291
	.value	0x2
	.long	.Ldebug_abbrev0
	.byte	0x8
	.uleb128 0x1
	.long	.LASF23
	.byte	0x1
	.long	.LASF24
	.long	.LASF25
	.quad	.Ltext0
	.quad	.Letext0
	.long	.Ldebug_line0
	.uleb128 0x2
	.byte	0x8
	.byte	0x7
	.long	.LASF0
	.uleb128 0x3
	.byte	0x4
	.byte	0x5
	.string	"int"
	.uleb128 0x2
	.byte	0x4
	.byte	0x7
	.long	.LASF1
	.uleb128 0x2
	.byte	0x8
	.byte	0x5
	.long	.LASF2
	.uleb128 0x2
	.byte	0x8
	.byte	0x5
	.long	.LASF3
	.uleb128 0x2
	.byte	0x1
	.byte	0x8
	.long	.LASF4
	.uleb128 0x2
	.byte	0x2
	.byte	0x7
	.long	.LASF5
	.uleb128 0x2
	.byte	0x1
	.byte	0x6
	.long	.LASF6
	.uleb128 0x2
	.byte	0x2
	.byte	0x5
	.long	.LASF7
	.uleb128 0x4
	.byte	0x8
	.byte	0x7
	.uleb128 0x5
	.byte	0x8
	.long	0x75
	.uleb128 0x2
	.byte	0x1
	.byte	0x6
	.long	.LASF8
	.uleb128 0x2
	.byte	0x8
	.byte	0x7
	.long	.LASF9
	.uleb128 0x6
	.long	.LASF11
	.byte	0x1
	.byte	0x6
	.long	0x8e
	.uleb128 0x2
	.byte	0x8
	.byte	0x4
	.long	.LASF10
	.uleb128 0x6
	.long	.LASF12
	.byte	0x1
	.byte	0x9
	.long	0x83
	.uleb128 0x6
	.long	.LASF13
	.byte	0x1
	.byte	0xd
	.long	0x42
	.uleb128 0x7
	.long	.LASF26
	.byte	0x4
	.byte	0x1
	.byte	0x1b
	.long	0xd4
	.uleb128 0x8
	.string	"row"
	.byte	0x1
	.byte	0x1c
	.long	0x65
	.byte	0x2
	.byte	0x23
	.uleb128 0x0
	.uleb128 0x9
	.long	.LASF14
	.byte	0x1
	.byte	0x1c
	.long	0x65
	.byte	0x2
	.byte	0x23
	.uleb128 0x2
	.byte	0x0
	.uleb128 0xa
	.byte	0x1
	.long	.LASF27
	.byte	0x1
	.byte	0x16
	.byte	0x1
	.long	0xfd
	.byte	0x3
	.long	0xfd
	.uleb128 0xb
	.string	"ptr"
	.byte	0x1
	.byte	0x16
	.long	0x6f
	.uleb128 0xc
	.long	.LASF15
	.byte	0x1
	.byte	0x17
	.long	0x103
	.byte	0x0
	.uleb128 0x5
	.byte	0x8
	.long	0x95
	.uleb128 0xd
	.long	0xa0
	.uleb128 0xe
	.byte	0x1
	.long	.LASF16
	.byte	0x1
	.byte	0xf
	.byte	0x1
	.long	0x6f
	.quad	.LFB12
	.quad	.LFE12
	.byte	0x2
	.byte	0x77
	.sleb128 8
	.long	0x145
	.uleb128 0xf
	.string	"ptr"
	.byte	0x1
	.byte	0xf
	.long	0x6f
	.byte	0x1
	.byte	0x55
	.uleb128 0xc
	.long	.LASF15
	.byte	0x1
	.byte	0x11
	.long	0xa0
	.byte	0x0
	.uleb128 0x10
	.long	0xd4
	.quad	.LFB13
	.quad	.LFE13
	.byte	0x2
	.byte	0x77
	.sleb128 8
	.long	0x16e
	.uleb128 0x11
	.long	0xe6
	.byte	0x1
	.byte	0x55
	.uleb128 0x12
	.long	0xf1
	.byte	0x0
	.uleb128 0x13
	.byte	0x1
	.long	.LASF17
	.byte	0x1
	.byte	0x20
	.byte	0x1
	.long	0x8e
	.quad	.LFB14
	.quad	.LFE14
	.long	.LLST2
	.long	0x25a
	.uleb128 0x14
	.string	"v1"
	.byte	0x1
	.byte	0x20
	.long	0x25a
	.long	.LLST3
	.uleb128 0x14
	.string	"v2"
	.byte	0x1
	.byte	0x20
	.long	0x25a
	.long	.LLST4
	.uleb128 0x15
	.long	.LASF18
	.byte	0x1
	.byte	0x20
	.long	0x8e
	.long	.LLST5
	.uleb128 0x14
	.string	"tsk"
	.byte	0x1
	.byte	0x20
	.long	0x6f
	.long	.LLST6
	.uleb128 0x16
	.string	"sum"
	.byte	0x1
	.byte	0x21
	.long	0x8e
	.long	.LLST7
	.uleb128 0xc
	.long	.LASF19
	.byte	0x1
	.byte	0x21
	.long	0x8e
	.uleb128 0xc
	.long	.LASF20
	.byte	0x1
	.byte	0x22
	.long	0x265
	.uleb128 0xc
	.long	.LASF21
	.byte	0x1
	.byte	0x23
	.long	0x6f
	.uleb128 0x17
	.long	.LASF22
	.byte	0x1
	.byte	0x25
	.long	0xfd
	.long	.LLST8
	.uleb128 0x18
	.long	0xd4
	.quad	.LBB5
	.long	.Ldebug_ranges0+0x0
	.byte	0x1
	.byte	0x25
	.long	0x234
	.uleb128 0x19
	.long	0x161
	.uleb128 0x1a
	.long	.Ldebug_ranges0+0x40
	.uleb128 0x12
	.long	0xf1
	.byte	0x0
	.byte	0x0
	.uleb128 0x1b
	.quad	.LBB11
	.quad	.LBE11
	.uleb128 0x1c
	.string	"i"
	.byte	0x1
	.byte	0x28
	.long	0x34
	.uleb128 0x1c
	.string	"tt"
	.byte	0x1
	.byte	0x31
	.long	0x274
	.byte	0x0
	.byte	0x0
	.uleb128 0x5
	.byte	0x8
	.long	0x260
	.uleb128 0xd
	.long	0x8e
	.uleb128 0x1d
	.long	0x8e
	.long	0x274
	.uleb128 0x1e
	.long	0x6c
	.byte	0x0
	.uleb128 0x5
	.byte	0x8
	.long	0x8e
	.uleb128 0x1f
	.string	"bs"
	.byte	0x1
	.byte	0x8
	.long	0x28f
	.byte	0x1
	.byte	0x9
	.byte	0x3
	.quad	bs
	.uleb128 0xd
	.long	0x34
	.byte	0x0
	.section	.debug_abbrev
	.uleb128 0x1
	.uleb128 0x11
	.byte	0x1
	.uleb128 0x25
	.uleb128 0xe
	.uleb128 0x13
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x1b
	.uleb128 0xe
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x10
	.uleb128 0x6
	.byte	0x0
	.byte	0x0
	.uleb128 0x2
	.uleb128 0x24
	.byte	0x0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0xe
	.byte	0x0
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x24
	.byte	0x0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0x8
	.byte	0x0
	.byte	0x0
	.uleb128 0x4
	.uleb128 0x24
	.byte	0x0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.byte	0x0
	.byte	0x0
	.uleb128 0x5
	.uleb128 0xf
	.byte	0x0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x6
	.uleb128 0x16
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x7
	.uleb128 0x13
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x8
	.uleb128 0xd
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x38
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x9
	.uleb128 0xd
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x38
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0xa
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x20
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0xb
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0xc
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0xd
	.uleb128 0x26
	.byte	0x0
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0xe
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x40
	.uleb128 0xa
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0xf
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x10
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x40
	.uleb128 0xa
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x11
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x12
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x31
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x13
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x40
	.uleb128 0x6
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x14
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x6
	.byte	0x0
	.byte	0x0
	.uleb128 0x15
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x6
	.byte	0x0
	.byte	0x0
	.uleb128 0x16
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x6
	.byte	0x0
	.byte	0x0
	.uleb128 0x17
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x6
	.byte	0x0
	.byte	0x0
	.uleb128 0x18
	.uleb128 0x1d
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x52
	.uleb128 0x1
	.uleb128 0x55
	.uleb128 0x6
	.uleb128 0x58
	.uleb128 0xb
	.uleb128 0x59
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x19
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x31
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1a
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x55
	.uleb128 0x6
	.byte	0x0
	.byte	0x0
	.uleb128 0x1b
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.byte	0x0
	.byte	0x0
	.uleb128 0x1c
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1d
	.uleb128 0x1
	.byte	0x1
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x1
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1e
	.uleb128 0x21
	.byte	0x0
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x1f
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.byte	0x0
	.section	.debug_pubnames,"",@progbits
	.long	0x61
	.value	0x2
	.long	.Ldebug_info0
	.long	0x295
	.long	0x108
	.string	"aligned_pointer"
	.long	0x145
	.string	"aligned_value_pointer"
	.long	0x16e
	.string	"asm_fast_sandwich_product"
	.long	0x27a
	.string	"bs"
	.long	0x0
	.section	.debug_aranges,"",@progbits
	.long	0x2c
	.value	0x2
	.long	.Ldebug_info0
	.byte	0x8
	.byte	0x0
	.value	0x0
	.value	0x0
	.quad	.Ltext0
	.quad	.Letext0-.Ltext0
	.quad	0x0
	.quad	0x0
	.section	.debug_ranges,"",@progbits
.Ldebug_ranges0:
	.quad	.LBB5-.Ltext0
	.quad	.LBE5-.Ltext0
	.quad	.LBB10-.Ltext0
	.quad	.LBE10-.Ltext0
	.quad	.LBB9-.Ltext0
	.quad	.LBE9-.Ltext0
	.quad	0x0
	.quad	0x0
	.quad	.LBB7-.Ltext0
	.quad	.LBE7-.Ltext0
	.quad	.LBB8-.Ltext0
	.quad	.LBE8-.Ltext0
	.quad	.LBB6-.Ltext0
	.quad	.LBE6-.Ltext0
	.quad	0x0
	.quad	0x0
	.section	.debug_str,"MS",@progbits,1
.LASF3:
	.string	"long long int"
.LASF20:
	.string	"vec1"
.LASF5:
	.string	"short unsigned int"
.LASF22:
	.string	"atsk"
.LASF14:
	.string	"column"
.LASF17:
	.string	"asm_fast_sandwich_product"
.LASF27:
	.string	"aligned_value_pointer"
.LASF4:
	.string	"unsigned char"
.LASF0:
	.string	"long unsigned int"
.LASF26:
	.string	"location_struct"
.LASF11:
	.string	"value_t"
.LASF10:
	.string	"double"
.LASF1:
	.string	"unsigned int"
.LASF24:
	.string	"packed_simulation_model.c"
.LASF9:
	.string	"long long unsigned int"
.LASF25:
	.string	"/home/axel/NewWeb/src"
.LASF18:
	.string	"factor"
.LASF15:
	.string	"mask"
.LASF23:
	.string	"GNU C 4.4.5"
.LASF13:
	.string	"size_T"
.LASF7:
	.string	"short int"
.LASF21:
	.string	"matrix_end"
.LASF12:
	.string	"aligned_value_t"
.LASF2:
	.string	"long int"
.LASF8:
	.string	"char"
.LASF6:
	.string	"signed char"
.LASF16:
	.string	"aligned_pointer"
.LASF19:
	.string	"last"
	.ident	"GCC: (Ubuntu/Linaro 4.4.4-14ubuntu5) 4.4.5"
	.section	.note.GNU-stack,"",@progbits
