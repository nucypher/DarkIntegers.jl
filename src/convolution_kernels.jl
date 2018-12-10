@inbounds function mul2_cyclic(x::NTuple{2, T}, y::NTuple{2, T}) where T
    a0, a1 = x
    b0, b1 = y

    z0 = a0 * b0 + a1 * b1
    z1 = (a0 + a1) * (b0 + b1) - z0

    (z0, z1)
end


@inbounds function mul2_negacyclic(x::NTuple{2, T}, y::NTuple{2, T}) where T
    a0, a1 = x
    b0, b1 = y

    t = a0 * (b0 + b1)
    z0 = t - (a0 + a1) * b1
    z1 = t + (a1 - a0) * b0

    (z0, z1)
end


@inbounds function mul4_cyclic_sb(x::NTuple{4, T}, y::NTuple{4, T}) where T
    a0, a1, a2, a3 = x
    b0, b1, b2, b3 = y

    z0 = a0*b0 + a1*b3 + a2*b2 + a3*b1
    z1 = a0*b1 + a1*b0 + a2*b3 + a3*b2
    z2 = a0*b2 + a1*b1 + a2*b0 + a3*b3
    z3 = a0*b3 + a1*b2 + a2*b1 + a3*b0

    (z0, z1, z2, z3)
end


@inbounds function mul4_negacyclic_sb(x::NTuple{4, T}, y::NTuple{4, T}) where T
    a0, a1, a2, a3 = x
    b0, b1, b2, b3 = y

    z0 = a0*b0 - a1*b3 - a2*b2 - a3*b1
    z1 = a0*b1 + a1*b0 - a2*b3 - a3*b2
    z2 = a0*b2 + a1*b1 + a2*b0 - a3*b3
    z3 = a0*b3 + a1*b2 + a2*b1 + a3*b0

    (z0, z1, z2, z3)
end


@inbounds function mul4_negacyclic_karatsuba(x::NTuple{4, T}, y::NTuple{4, T}) where T
    a0, a1, a2, a3 = x
    b0, b1, b2, b3 = y

    # Unrolled 2-stage Karatsuba
    t0 = a0*b0
    t4 = a1*b1
    t5 = t0 + t4 + (a1 - a0)*(b0 - b1)
    t6 = a2 - a0
    t7 = a3 - a1
    t8 = b0 - b2
    t9 = b1 - b3
    t10 = t6*t8
    t14 = t7*t9
    t16 = a2*b2
    t20 = a3*b3
    t21 = t16 + t20 + (a3 - a2)*(b2 - b3)
    t23 = t10 + t14 + t21 + t5 + (t7 - t6)*(t8 - t9)
    z0 = t0 - t14 - t16 - t20 - t4
    z1 = t5 - t21
    z2 = t0 + t10 + t16 + t4 - t20
    z3 = t23

    (z0, z1, z2, z3)
end


@inbounds function mul4_cyclic_karatsuba(x::NTuple{4, T}, y::NTuple{4, T}) where T
    a0, a1, a2, a3 = x
    b0, b1, b2, b3 = y

    # Unrolled 2-stage Karatsuba
    t0 = a0*b0
    t4 = a1*b1
    t5 = t0 + t4 + (a1 - a0)*(b0 - b1)
    t6 = a2 - a0
    t7 = a3 - a1
    t8 = b0 - b2
    t9 = b1 - b3
    t10 = t6*t8
    t14 = t7*t9
    t16 = a2*b2
    t20 = a3*b3
    t21 = t16 + t20 + (a3 - a2)*(b2 - b3)
    t23 = t10 + t14 + t21 + t5 + (t7 - t6)*(t8 - t9)
    z0 = t0 + t14 + t16 + t20 + t4
    z1 = t21 + t5
    z2 = t0 + t10 + t16 + t20 + t4
    z3 = t23

    (z0, z1, z2, z3)
end


@inbounds function mul8_negacyclic_karatsuba(x::NTuple{8, T}, y::NTuple{8, T}) where T

    a0, a1, a2, a3, a4, a5, a6, a7 = x
    b0, b1, b2, b3, b4, b5, b6, b7 = y

    # Unrolled 3-stage Karatsuba
    t0 = a0*b0
    t4 = a1*b1
    t5 = t0 + t4 + (a1 - a0)*(b0 - b1)
    t6 = a2 - a0
    t7 = a3 - a1
    t8 = b0 - b2
    t9 = b1 - b3
    t10 = t6*t8
    t14 = t7*t9
    t16 = a2*b2
    t20 = a3*b3
    t21 = t16 + t20 + (a3 - a2)*(b2 - b3)
    t23 = t10 + t14 + t21 + t5 + (t7 - t6)*(t8 - t9)
    t25 = t0 + t10 + t16 + t4
    t26 = t14 + t16 + t20 + t4
    t27 = a4 - a0
    t28 = a5 - a1
    t29 = a6 - a2
    t30 = a7 - a3
    t31 = b0 - b4
    t32 = b1 - b5
    t33 = b2 - b6
    t34 = b3 - b7
    t35 = t27*t31
    t39 = t28*t32
    t40 = t35 + t39 + (t28 - t27)*(t31 - t32)
    t41 = t29 - t27
    t42 = t30 - t28
    t43 = t31 - t33
    t44 = t32 - t34
    t45 = t41*t43
    t49 = t42*t44
    t51 = t29*t33
    t55 = t30*t34
    t56 = t51 + t55 + (t30 - t29)*(t33 - t34)
    t62 = a4*b4
    t66 = a5*b5
    t67 = t62 + t66 + (a5 - a4)*(b4 - b5)
    t68 = a6 - a4
    t69 = a7 - a5
    t70 = b4 - b6
    t71 = b5 - b7
    t72 = t68*t70
    t76 = t69*t71
    t78 = a6*b6
    t82 = a7*b7
    t83 = t78 + t82 + (a7 - a6)*(b6 - b7)
    t85 = t67 + t72 + t76 + t83 + (t69 - t68)*(t70 - t71)
    t87 = t62 + t66 + t72 + t78
    t88 = t66 + t76 + t78 + t82
    t92 = t23 + t40 + t45 + t49 + t56 + t85 + (t42 - t41)*(t43 - t44)

    z0 = t0 - t26 - t39 - t49 - t51 - t55 - t62 - t88
    z1 = t5 - t21 - t56 - t67 - t83
    z2 = t25 - t20 - t55 - t82 - t87
    z3 = t23 - t85
    z4 = t0 + t26 + t35 + t62 - t88
    z5 = t21 + t40 + t5 + t67 - t83
    z6 = t20 + t25 + t35 + t39 + t45 + t51 + t87 - t82
    z7 = t92

    (z0, z1, z2, z3, z4, z5, z6, z7)
end


@inbounds function mul8_cyclic_karatsuba(x::NTuple{8, T}, y::NTuple{8, T}) where T

    a0, a1, a2, a3, a4, a5, a6, a7 = x
    b0, b1, b2, b3, b4, b5, b6, b7 = y

    # Unrolled 3-stage Karatsuba
    t0 = a0*b0
    t4 = a1*b1
    t5 = t0 + t4 + (a1 - a0)*(b0 - b1)
    t6 = a2 - a0
    t7 = a3 - a1
    t8 = b0 - b2
    t9 = b1 - b3
    t10 = t6*t8
    t14 = t7*t9
    t16 = a2*b2
    t20 = a3*b3
    t21 = t16 + t20 + (a3 - a2)*(b2 - b3)
    t23 = t10 + t14 + t21 + t5 + (t7 - t6)*(t8 - t9)
    t25 = t0 + t10 + t16 + t4
    t26 = t14 + t16 + t20 + t4
    t27 = a4 - a0
    t28 = a5 - a1
    t29 = a6 - a2
    t30 = a7 - a3
    t31 = b0 - b4
    t32 = b1 - b5
    t33 = b2 - b6
    t34 = b3 - b7
    t35 = t27*t31
    t39 = t28*t32
    t40 = t35 + t39 + (t28 - t27)*(t31 - t32)
    t41 = t29 - t27
    t42 = t30 - t28
    t43 = t31 - t33
    t44 = t32 - t34
    t45 = t41*t43
    t49 = t42*t44
    t51 = t29*t33
    t55 = t30*t34
    t56 = t51 + t55 + (t30 - t29)*(t33 - t34)
    t62 = a4*b4
    t66 = a5*b5
    t67 = t62 + t66 + (a5 - a4)*(b4 - b5)
    t68 = a6 - a4
    t69 = a7 - a5
    t70 = b4 - b6
    t71 = b5 - b7
    t72 = t68*t70
    t76 = t69*t71
    t78 = a6*b6
    t82 = a7*b7
    t83 = t78 + t82 + (a7 - a6)*(b6 - b7)
    t85 = t67 + t72 + t76 + t83 + (t69 - t68)*(t70 - t71)
    t87 = t62 + t66 + t72 + t78
    t88 = t66 + t76 + t78 + t82
    t92 = t23 + t40 + t45 + t49 + t56 + t85 + (t42 - t41)*(t43 - t44)
    z0 = t0 + t26 + t39 + t49 + t51 + t55 + t62 + t88
    z1 = t21 + t5 + t56 + t67 + t83
    z2 = t20 + t25 + t55 + t82 + t87
    z3 = t23 + t85
    z4 = t0 + t26 + t35 + t62 + t88
    z5 = t21 + t40 + t5 + t67 + t83
    z6 = t20 + t25 + t35 + t39 + t45 + t51 + t82 + t87
    z7 = t92

    (z0, z1, z2, z3, z4, z5, z6, z7)
end
