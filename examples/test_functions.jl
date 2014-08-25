include("nsga2.jl")

# benchmarking for nsga2


# def zdt1(individual):
#     """ZDT1 multiobjective function.
#     
#     :math:`g(\\mathbf{x}) = 1 + \\frac{9}{n-1}\\sum_{i=2}^n x_i`
#     
#     :math:`f_{\\text{ZDT1}1}(\\mathbf{x}) = x_1`
#     
#     :math:`f_{\\text{ZDT1}2}(\\mathbf{x}) = g(\\mathbf{x})\\left[1 - \\sqrt{\\frac{x_1}{g(\\mathbf{x})}}\\right]`
#     """
#     g  = 1.0 + 9.0*sum(individual[1:])/(len(individual)-1)
#     f1 = individual[0]
#     f2 = g * (1 - sqrt(f1/g))
#     return f1, f2

function zdt1(individual::Individual)
  # ZDT1 multiobjective function
  g::FloatingPoint = 1.0 + (9.0 * sum(individual.genes[2:end])/(length(individual.genes) -1))
  f1::FloatingPoint = individual.genes[1]
  f2::FloatingPoint = g * (1 - sqrt(f1 / g))
  (f1, f2)
end


# def zdt2(individual):
#     """ZDT2 multiobjective function.
#     
#     :math:`g(\\mathbf{x}) = 1 + \\frac{9}{n-1}\\sum_{i=2}^n x_i`
#     
#     :math:`f_{\\text{ZDT2}1}(\\mathbf{x}) = x_1`
#     
#     :math:`f_{\\text{ZDT2}2}(\\mathbf{x}) = g(\\mathbf{x})\\left[1 - \\left(\\frac{x_1}{g(\\mathbf{x})}\\right)^2\\right]`
#     
#     """
# 
#     g  = 1.0 + 9.0*sum(individual[1:])/(len(individual)-1)
#     f1 = individual[0]
#     f2 = g * (1 - (f1/g)**2)
#     return f1, f2

function zdt2(individual::Individual)
  # ZDT2 multiobjective function
  g::FloatingPoint = 1.0 + (9.0 * sum(individual.genes[2:end])/(length(individual.genes) -1))
  f1::FloatingPoint = individual.genes[1]
  f2::FloatingPoint = g * (1 - (f1 / g) ^ 2)  # verify the exponentiation in julia
  (f1, f2)
end


# def zdt3(individual):
#     """ZDT3 multiobjective function.
# 
#     :math:`g(\\mathbf{x}) = 1 + \\frac{9}{n-1}\\sum_{i=2}^n x_i`
# 
#     :math:`f_{\\text{ZDT3}1}(\\mathbf{x}) = x_1`
# 
#     :math:`f_{\\text{ZDT3}2}(\\mathbf{x}) = g(\\mathbf{x})\\left[1 - \\sqrt{\\frac{x_1}{g(\\mathbf{x})}} - \\frac{x_1}{g(\\mathbf{x})}\\sin(10\\pi x_1)\\right]`
# 
#     """
# 
#     g  = 1.0 + 9.0*sum(individual[1:])/(len(individual)-1)
#     f1 = individual[0]
#     f2 = g * (1 - sqrt(f1/g) - f1/g * sin(10*pi*f1))
#     return f1, f2

function zdt3(individual::Individual)
  g::FloatingPoint = 1.0 + (9.0 * sum(individual.genes[2:end])/(length(individual.genes) -1))
  f1::FloatingPoint = individual.genes[1]
  f2::FloatingPoint = g * (1 - (f1 / g) ^ 2)  # verify the exponentiation in julia
  (f1, f2)
end


# def zdt4(individual):
#     """ZDT4 multiobjective function.
#     
#     :math:`g(\\mathbf{x}) = 1 + 10(n-1) + \\sum_{i=2}^n \\left[ x_i^2 - 10\\cos(4\\pi x_i) \\right]`
# 
#     :math:`f_{\\text{ZDT4}1}(\\mathbf{x}) = x_1`
#     
#     :math:`f_{\\text{ZDT4}2}(\\mathbf{x}) = g(\\mathbf{x})\\left[ 1 - \\sqrt{x_1/g(\\mathbf{x})} \\right]`
#     
#     """
#     g  = 1 + 10*(len(individual)-1) + sum(xi**2 - 10*cos(4*pi*xi) for xi in individual[1:])
#     f1 = individual[0]
#     f2 = g * (1 - sqrt(f1/g))
#     return f1, f2

function zdt4(individual::Individual)
  # ZDT4 multiobjective function
  g::FloatingPoint = 1.0 + (10.0 * (length(individual.genes) - 1) + sum()
  f1::FloatingPoint = individual.genes[1]
  f2::FloatingPoint = g * (1 - (f1 / g) ^ 2)  # verify the exponentiation in julia
  (f1, f2)
end


# def zdt6(individual):
#     """ZDT6 multiobjective function.
#     
#     :math:`g(\\mathbf{x}) = 1 + 9 \\left[ \\left(\\sum_{i=2}^n x_i\\right)/(n-1) \\right]^{0.25}`
#     
#     :math:`f_{\\text{ZDT6}1}(\\mathbf{x}) = 1 - \\exp(-4x_1)\\sin^6(6\\pi x_1)`
#     
#     :math:`f_{\\text{ZDT6}2}(\\mathbf{x}) = g(\\mathbf{x}) \left[ 1 - (f_{\\text{ZDT6}1}(\\mathbf{x})/g(\\mathbf{x}))^2 \\right]`
#     
#     """
#     g  = 1 + 9 * (sum(individual[1:]) / (len(individual)-1))**0.25
#     f1 = 1 - exp(-4*individual[0]) * sin(6*pi*individual[0])**6
#     f2 = g * (1 - (f1/g)**2)
#     return f1, f2

function zdt6(individual::Individual)
  # ZDT6 multiobjective function
  g::FloatingPoint = 1.0 + (9.0 * sum(individual.genes[2:end])/((length(individual.genes) -1))^0.25)
  f1::FloatingPoint = 1 - (exp(-4 * individual.genes[1])) * (sin(6 * pi * individual[1])^6))
  f2::FloatingPoint = g * (1 - (f1 / g) ^ 2)  # verify the exponentiation in julia
  (f1, f2)
end


# def dtlz1(individual, obj):
#     """DTLZ1 mutliobjective function. It returns a tuple of *obj* values. 
#     The individual must have at least *obj* elements.
#     From: K. Deb, L. Thiele, M. Laumanns and E. Zitzler. Scalable Multi-Objective 
#     Optimization Test Problems. CEC 2002, p. 825 - 830, IEEE Press, 2002.
# 
#     :math:`g(\\mathbf{x}_m) = 100\\left(|\\mathbf{x}_m| + \sum_{x_i \in \\mathbf{x}_m}\\left((x_i - 0.5)^2 - \\cos(20\pi(x_i - 0.5))\\right)\\right)`
# 
#     :math:`f_{\\text{DTLZ1}1}(\\mathbf{x}) = \\frac{1}{2} (1 + g(\\mathbf{x}_m)) \\prod_{i=1}^{m-1}x_i`
#     
#     :math:`f_{\\text{DTLZ1}2}(\\mathbf{x}) = \\frac{1}{2} (1 + g(\\mathbf{x}_m)) (1-x_{m-1}) \\prod_{i=1}^{m-2}x_i`
#     
#     :math:`\\ldots`
#     
#     :math:`f_{\\text{DTLZ1}m-1}(\\mathbf{x}) = \\frac{1}{2} (1 + g(\\mathbf{x}_m)) (1 - x_2) x_1`
#     
#     :math:`f_{\\text{DTLZ1}m}(\\mathbf{x}) = \\frac{1}{2} (1 - x_1)(1 + g(\\mathbf{x}_m))`
#     
#     Where :math:`m` is the number of objectives and :math:`\\mathbf{x}_m` is a
#     vector of the remaining attributes :math:`[x_m~\\ldots~x_n]` of the
#     individual in :math:`n > m` dimensions.
#     
#     """
#     g = 100 * (len(individual[obj-1:]) + sum((xi-0.5)**2 - cos(20*pi*(xi-0.5)) for xi in individual[obj-1:]))
#     f = [0.5 * reduce(mul, individual[:obj-1], 1) * (1 + g)]
#     f.extend(0.5 * reduce(mul, individual[:m], 1) * (1 - individual[m]) * (1 + g) for m in reversed(xrange(obj-1)))
#     return f
# 
# def dtlz2(individual, obj):
#     """DTLZ2 mutliobjective function. It returns a tuple of *obj* values. 
#     The individual must have at least *obj* elements.
#     From: K. Deb, L. Thiele, M. Laumanns and E. Zitzler. Scalable Multi-Objective 
#     Optimization Test Problems. CEC 2002, p. 825 - 830, IEEE Press, 2002.
#     
#     :math:`g(\\mathbf{x}_m) = \\sum_{x_i \in \\mathbf{x}_m} (x_i - 0.5)^2`
#     
#     :math:`f_{\\text{DTLZ2}1}(\\mathbf{x}) = (1 + g(\\mathbf{x}_m)) \\prod_{i=1}^{m-1} \\cos(0.5x_i\pi)`
#     
#     :math:`f_{\\text{DTLZ2}2}(\\mathbf{x}) = (1 + g(\\mathbf{x}_m)) \\sin(0.5x_{m-1}\pi ) \\prod_{i=1}^{m-2} \\cos(0.5x_i\pi)`
#     
#     :math:`\\ldots`
#     
#     :math:`f_{\\text{DTLZ2}m}(\\mathbf{x}) = (1 + g(\\mathbf{x}_m)) \\sin(0.5x_{1}\pi )`
#     
#     Where :math:`m` is the number of objectives and :math:`\\mathbf{x}_m` is a
#     vector of the remaining attributes :math:`[x_m~\\ldots~x_n]` of the
#     individual in :math:`n > m` dimensions.
#     """
#     xc = individual[:obj-1]
#     xm = individual[obj-1:]
#     g = sum((xi-0.5)**2 for xi in xm)
#     f = [(1.0+g) *  reduce(mul, (cos(0.5*xi*pi) for xi in xc), 1.0)]
#     f.extend((1.0+g) * reduce(mul, (cos(0.5*xi*pi) for xi in xc[:m]), 1) * sin(0.5*xc[m]*pi) for m in range(obj-2, -1, -1))
# 
#     return f
# 
# def dtlz3(individual, obj):
#     """DTLZ3 mutliobjective function. It returns a tuple of *obj* values. 
#     The individual must have at least *obj* elements.
#     From: K. Deb, L. Thiele, M. Laumanns and E. Zitzler. Scalable Multi-Objective 
#     Optimization Test Problems. CEC 2002, p. 825 - 830, IEEE Press, 2002.
#     
#     :math:`g(\\mathbf{x}_m) = 100\\left(|\\mathbf{x}_m| + \sum_{x_i \in \\mathbf{x}_m}\\left((x_i - 0.5)^2 - \\cos(20\pi(x_i - 0.5))\\right)\\right)`
#     
#     :math:`f_{\\text{DTLZ3}1}(\\mathbf{x}) = (1 + g(\\mathbf{x}_m)) \\prod_{i=1}^{m-1} \\cos(0.5x_i\pi)`
#     
#     :math:`f_{\\text{DTLZ3}2}(\\mathbf{x}) = (1 + g(\\mathbf{x}_m)) \\sin(0.5x_{m-1}\pi ) \\prod_{i=1}^{m-2} \\cos(0.5x_i\pi)`
#     
#     :math:`\\ldots`
#     
#     :math:`f_{\\text{DTLZ3}m}(\\mathbf{x}) = (1 + g(\\mathbf{x}_m)) \\sin(0.5x_{1}\pi )`
#     
#     Where :math:`m` is the number of objectives and :math:`\\mathbf{x}_m` is a
#     vector of the remaining attributes :math:`[x_m~\\ldots~x_n]` of the
#     individual in :math:`n > m` dimensions.
#     """
#     xc = individual[:obj-1]
#     xm = individual[obj-1:]
#     g = 100 * (len(xm) + sum((xi-0.5)**2 - cos(20*pi*(xi-0.5)) for xi in xm))
#     f = [(1.0+g) *  reduce(mul, (cos(0.5*xi*pi) for xi in xc), 1.0)]
#     f.extend((1.0+g) * reduce(mul, (cos(0.5*xi*pi) for xi in xc[:m]), 1) * sin(0.5*xc[m]*pi) for m in range(obj-2, -1, -1))
#     return f
# 
# def dtlz4(individual, obj, alpha):
#     """DTLZ4 mutliobjective function. It returns a tuple of *obj* values. The
#     individual must have at least *obj* elements. The *alpha* parameter allows
#     for a meta-variable mapping in :func:`dtlz2` :math:`x_i \\rightarrow
#     x_i^\\alpha`, the authors suggest :math:`\\alpha = 100`.
#     From: K. Deb, L. Thiele, M. Laumanns and E. Zitzler. Scalable Multi-Objective 
#     Optimization Test Problems. CEC 2002, p. 825 - 830, IEEE Press, 2002.
#     
#     :math:`g(\\mathbf{x}_m) = \\sum_{x_i \in \\mathbf{x}_m} (x_i - 0.5)^2`
#     
#     :math:`f_{\\text{DTLZ4}1}(\\mathbf{x}) = (1 + g(\\mathbf{x}_m)) \\prod_{i=1}^{m-1} \\cos(0.5x_i^\\alpha\pi)`
#     
#     :math:`f_{\\text{DTLZ4}2}(\\mathbf{x}) = (1 + g(\\mathbf{x}_m)) \\sin(0.5x_{m-1}^\\alpha\pi ) \\prod_{i=1}^{m-2} \\cos(0.5x_i^\\alpha\pi)`
#     
#     :math:`\\ldots`
#     
#     :math:`f_{\\text{DTLZ4}m}(\\mathbf{x}) = (1 + g(\\mathbf{x}_m)) \\sin(0.5x_{1}^\\alpha\pi )`
#     
#     Where :math:`m` is the number of objectives and :math:`\\mathbf{x}_m` is a
#     vector of the remaining attributes :math:`[x_m~\\ldots~x_n]` of the
#     individual in :math:`n > m` dimensions.
#     """
#     xc = individual[:obj-1]
#     xm = individual[obj-1:]
#     g = sum((xi-0.5)**2 for xi in xm)
#     f = [(1.0+g) *  reduce(mul, (cos(0.5*xi**alpha*pi) for xi in xc), 1.0)]
#     f.extend((1.0+g) * reduce(mul, (cos(0.5*xi**alpha*pi) for xi in xc[:m]), 1) * sin(0.5*xc[m]**alpha*pi) for m in range(obj-2, -1, -1))
#     return f
# 
# def fonseca(individual):
#     """Fonseca and Fleming's multiobjective function.
#     From: C. M. Fonseca and P. J. Fleming, "Multiobjective optimization and
#     multiple constraint handling with evolutionary algorithms -- Part II:
#     Application example", IEEE Transactions on Systems, Man and Cybernetics,
#     1998.
#     
#     :math:`f_{\\text{Fonseca}1}(\\mathbf{x}) = 1 - e^{-\\sum_{i=1}^{3}(x_i - \\frac{1}{\\sqrt{3}})^2}`
#     
#     :math:`f_{\\text{Fonseca}2}(\\mathbf{x}) = 1 - e^{-\\sum_{i=1}^{3}(x_i + \\frac{1}{\\sqrt{3}})^2}`
#     """
#     f_1 = 1 - exp(-sum((xi - 1/sqrt(3))**2 for xi in individual[:3]))
#     f_2 = 1 - exp(-sum((xi + 1/sqrt(3))**2 for xi in individual[:3]))
#     return f_1, f_2
# 
# def poloni(individual):
#     """Poloni's multiobjective function on a two attribute *individual*. From:
#     C. Poloni, "Hybrid GA for multi objective aerodynamic shape optimization",
#     in Genetic Algorithms in Engineering and Computer Science, 1997.
#     
#     :math:`A_1 = 0.5 \\sin (1) - 2 \\cos (1) + \\sin (2) - 1.5 \\cos (2)`
# 
#     :math:`A_2 = 1.5 \\sin (1) - \\cos (1) + 2 \\sin (2) - 0.5 \\cos (2)`
# 
#     :math:`B_1 = 0.5 \\sin (x_1) - 2 \\cos (x_1) + \\sin (x_2) - 1.5 \\cos (x_2)`
# 
#     :math:`B_2 = 1.5 \\sin (x_1) - cos(x_1) + 2 \\sin (x_2) - 0.5 \\cos (x_2)`
#     
#     :math:`f_{\\text{Poloni}1}(\\mathbf{x}) = 1 + (A_1 - B_1)^2 + (A_2 - B_2)^2`
#     
#     :math:`f_{\\text{Poloni}2}(\\mathbf{x}) = (x_1 + 3)^2 + (x_2 + 1)^2`
#     """
#     x_1 = individual[0]
#     x_2 = individual[1]
#     A_1 = 0.5 * sin(1) - 2 * cos(1) + sin(2) - 1.5 * cos(2)
#     A_2 = 1.5 * sin(1) - cos(1) + 2 * sin(2) - 0.5 * cos(2)
#     B_1 = 0.5 * sin(x_1) - 2 * cos(x_1) + sin(x_2) - 1.5 * cos(x_2)
#     B_2 = 1.5 * sin(x_1) - cos(x_1) + 2 * sin(x_2) - 0.5 * cos(x_2)
#     return 1 + (A_1 - B_1)**2 + (A_2 - B_2)**2, (x_1 + 3)**2 + (x_2 + 1)**2
# 
