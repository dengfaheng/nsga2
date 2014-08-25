include("nsga2.jl")

#benchmarking for nsga2


#def zdt1(individual):
#    ZDT1 multiobjective function
#    g  = 1.0 + 9.0*sum(individual[1:])/(len(individual)-1)
#    f1 = individual[0]
#    f2 = g * (1 - sqrt(f1/g))
#    return f1, f2

function zdt1(individual::Individual)
  #ZDT1 multiobjective function
  g::FloatingPoint = 1.0 + (9.0 * sum(individual.genes[2:end])/(length(individual.genes) -1))
  f1::FloatingPoint = individual.genes[1]
  f2::FloatingPoint = g * (1 - sqrt(f1 / g))
  (f1, f2)
end


#def zdt2(individual):
#    ZDT2 multiobjective function
#    g  = 1.0 + 9.0*sum(individual[1:])/(len(individual)-1)
#    f1 = individual[0]
#    f2 = g * (1 - (f1/g)**2)
#    return f1, f2

function zdt2(individual::Individual)
  #ZDT2 multiobjective function
  g::FloatingPoint = 1.0 + (9.0 * sum(individual.genes[2:end])/(length(individual.genes) -1))
  f1::FloatingPoint = individual.genes[1]
  f2::FloatingPoint = g * (1 - (f1 / g) ^ 2)
  (f1, f2)
end


#def zdt3(individual):
#    ZDT3 multiobjective function
#    g  = 1.0 + 9.0*sum(individual[1:])/(len(individual)-1)
#    f1 = individual[0]
#    f2 = g * (1 - sqrt(f1/g) - f1/g * sin(10*pi*f1))
#    return f1, f2

function zdt3(individual::Individual)
  g::FloatingPoint = 1.0 + (9.0 * sum(individual.genes[2:end])/(length(individual.genes) -1))
  f1::FloatingPoint = individual.genes[1]
  f2::FloatingPoint = g * (1 - sqrt(f1 / g) - (f1/g * sin(10 * pi * f1)))
  (f1, f2)
end


#def zdt4(individual):
#    ZDT4 multiobjective function
#    g  = 1 + 10*(len(individual)-1) + sum(xi**2 - 10*cos(4*pi*xi) for xi in individual[1:])
#    f1 = individual[0]
#    f2 = g * (1 - sqrt(f1/g))
#    return f1, f2

function zdt4(individual::Individual)
  # ZDT4 multiobjective function
  #    g  = 1 + 10*(len(individual)-1) + sum(xi**2 - 10*cos(4*pi*xi) for xi in individual[1:])
  g::FloatingPoint = 1.0 + (10.0 * (length(individual.genes) - 1)) + sum(map(x -> x^2.0 - 10 *cos(4 * pi * x), individual.genes[2:end]))
  f1::FloatingPoint = individual.genes[1]
  f2::FloatingPoint = g * (1 - sqrt(f1 / g))
  (f1, f2)
end


#def zdt6(individual):
#    ZDT6 multiobjective function
#    g  = 1 + 9 * (sum(individual[1:]) / (len(individual)-1))**0.25
#    f1 = 1 - exp(-4*individual[0]) * sin(6*pi*individual[0])**6
#    f2 = g * (1 - (f1/g)**2)
#    return f1, f2

function zdt6(individual::Individual)
  #ZDT6 multiobjective function
  g::FloatingPoint = 1.0 + (9.0 * sum(individual.genes[2:end])/((length(individual.genes) -1)) ^ 0.25)
  f1::FloatingPoint = 1 - (exp(-4 * individual.genes[1])) * (sin(6 * pi * individual[1]) ^ 6.0))
  f2::FloatingPoint = g * (1 - (f1 / g) ^ 2.0)
  (f1, f2)
end


#def dtlz1(individual, obj):
#    DTLZ1 mutliobjective function. It returns a tuple of *obj* values
#    The individual must have at least *obj* elements.
#    From: K. Deb, L. Thiele, M. Laumanns and E. Zitzler. Scalable Multi-Objective
#    Optimization Test Problems. CEC 2002, p. 825 - 830, IEEE Press, 2002.
#    g = 100 * (len(individual[obj-1:]) + sum((xi-0.5)**2 - cos(20*pi*(xi-0.5)) for xi in individual[obj-1:]))
#    f = [0.5 * reduce(mul, individual[:obj-1], 1) * (1 + g)]
#    f.extend(0.5 * reduce(mul, individual[:m], 1) * (1 - individual[m]) * (1 + g) for m in reversed(xrange(obj-1)))
#    return f

function dtlz1(individual, obj)

end


#def dtlz2(individual, obj):
#    DTLZ2 mutliobjective function. It returns a tuple of *obj* values.
#    The individual must have at least *obj* elements.
#    From: K. Deb, L. Thiele, M. Laumanns and E. Zitzler. Scalable Multi-Objective
#    Optimization Test Problems. CEC 2002, p. 825 - 830, IEEE Press, 2002.
#    xc = individual[:obj-1]
#    xm = individual[obj-1:]
#    g = sum((xi-0.5)**2 for xi in xm)
#    f = [(1.0+g) *  reduce(mul, (cos(0.5*xi*pi) for xi in xc), 1.0)]
#    f.extend((1.0+g) * reduce(mul, (cos(0.5*xi*pi) for xi in xc[:m]), 1) * sin(0.5*xc[m]*pi) for m in range(obj-2, -1, -1))
#
#    return f



#def dtlz3(individual, obj):
#    DTLZ3 mutliobjective function. It returns a tuple of *obj* values.
#    The individual must have at least *obj* elements.
#    From: K. Deb, L. Thiele, M. Laumanns and E. Zitzler. Scalable Multi-Objective
#    Optimization Test Problems. CEC 2002, p. 825 - 830, IEEE Press, 2002.
#    xc = individual[:obj-1]
#    xm = individual[obj-1:]
#    g = 100 * (len(xm) + sum((xi-0.5)**2 - cos(20*pi*(xi-0.5)) for xi in xm))
#    f = [(1.0+g) *  reduce(mul, (cos(0.5*xi*pi) for xi in xc), 1.0)]
#    f.extend((1.0+g) * reduce(mul, (cos(0.5*xi*pi) for xi in xc[:m]), 1) * sin(0.5*xc[m]*pi) for m in range(obj-2, -1, -1))
#    return f



#def dtlz4(individual, obj, alpha):
#    DTLZ4 mutliobjective function. It returns a tuple of *obj* values. The
#    individual must have at least *obj* elements. The *alpha* parameter allows
#    for a meta-variable mapping in :func:`dtlz2` :math:`x_i \\rightarrow
#    x_i^\\alpha`, the authors suggest :math:`\\alpha = 100`.
#    From: K. Deb, L. Thiele, M. Laumanns and E. Zitzler. Scalable Multi-Objective
#    Optimization Test Problems. CEC 2002, p. 825 - 830, IEEE Press, 2002.
#    xc = individual[:obj-1]
#    xm = individual[obj-1:]
#    g = sum((xi-0.5)**2 for xi in xm)
#    f = [(1.0+g) *  reduce(mul, (cos(0.5*xi**alpha*pi) for xi in xc), 1.0)]
#    f.extend((1.0+g) * reduce(mul, (cos(0.5*xi**alpha*pi) for xi in xc[:m]), 1) * sin(0.5*xc[m]**alpha*pi) for m in range(obj-2, -1, -1))
#    return f

function dtlz4(individual, obj, alpha)


end

#def fonseca(individual):
#    Fonseca and Fleming's multiobjective function.
#    From: C. M. Fonseca and P. J. Fleming, "Multiobjective optimization and
#    multiple constraint handling with evolutionary algorithms -- Part II:
#    Application example", IEEE Transactions on Systems, Man and Cybernetics,
#    1998.
#    f_1 = 1 - exp(-sum((xi - 1/sqrt(3))**2 for xi in individual[:3]))
#    f_2 = 1 - exp(-sum((xi + 1/sqrt(3))**2 for xi in individual[:3]))
#    return f_1, f_2

function fonseca(individual::Individual)
  #    Fonseca and Fleming's multiobjective function.
  #    From: C. M. Fonseca and P. J. Fleming, "Multiobjective optimization and
  # multiple constraint handling with evolutionary algorithms -- Part II:
  # Application example", IEEE Transactions on Systems, Man and Cybernetics,
  # 1998.
  f_1 = 1 - exp(-1 * sum(map(x-> (x - 1/sqrt(3)) ^ 2, individual.genes[1:3])))
  f_1 = 1 - exp(-1 * sum(map(x-> (x + 1/sqrt(3)) ^ 2, individual.genes[1:3])))
  (f1, f2)
end


#def poloni(individual):
#    """Poloni's multiobjective function on a two attribute *individual*. From:
#    C. Poloni, "Hybrid GA for multi objective aerodynamic shape optimization",
#    in Genetic Algorithms in Engineering and Computer Science, 1997.
#    """
#    x_1 = individual[0]
#    x_2 = individual[1]
#    A_1 = 0.5 * sin(1) - 2 * cos(1) + sin(2) - 1.5 * cos(2)
#    A_2 = 1.5 * sin(1) - cos(1) + 2 * sin(2) - 0.5 * cos(2)
#    B_1 = 0.5 * sin(x_1) - 2 * cos(x_1) + sin(x_2) - 1.5 * cos(x_2)
#    B_2 = 1.5 * sin(x_1) - cos(x_1) + 2 * sin(x_2) - 0.5 * cos(x_2)
#    return 1 + (A_1 - B_1)**2 + (A_2 - B_2)**2, (x_1 + 3)**2 + (x_2 + 1)**2

function poloni(individual::Individual)
  # Poloni's multiobjective function on a two attribute *individual*. From:
  # C. Poloni, "Hybrid GA for multi objective aerodynamic shape optimization",
  # in Genetic Algorithms in Engineering and Computer Science, 1997.
  x_1 = individual.genes[1]
  x_2 = individual.genes[2]
  A_1 = 0.5 * sin(1) - 2 * cos(1) + sin(2) - 1.5 * cos(2)
  A_2 = 1.5 * sin(1) - cos(1) + 2 * sin(2) - 0.5 * cos(2)
  B_1 = 0.5 * sin(x_1) - 2 * cos(x_1) + sin(x_2) - 1.5 * cos(x_2)
  B_2 = 1.5 * sin(x_1) - cos(x_1) + 2 * sin(x_2) - 0.5 * cos(x_2)
  1 + (A_1 - B_1) ^ 2 + (A_2 - B_2) ^ 2, (x_1 + 3) ^ 2 + (x_2 + 1) ^ 2
end

