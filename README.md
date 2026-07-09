Maybe add info about read and write ports count? That's why we prefer v_st over v_stu_p3 when possible.

> subnormals: Zen handles them at ~full speed while Intel takes microcode assists in several cases - largely moot for us because we flush denormals to zero (see target hardware)
Are you sure?


Our platforms have different calling conventions for vector registers, we may say about it. On linux all xmm are scratch, on windows some registers should not be changed by function. Parameters passing in registers is also different.

