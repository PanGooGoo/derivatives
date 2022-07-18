from scipy import integrate
from numpy import inf
import numpy

'''
@note: AntonovKonikovSpector (with Rho=0).
Particular care should be taken about the ATM strike.
This code should be further reviewed as bugs might still be present.
'''

class AKSZeroCorrelationSABR():
    
    def __init__(self, alpha, beta, nu, expiry, strike, fwd, riemann = False, debug = False):
        self.alpha = alpha
        self.beta = beta
        self.nu = nu
        self.expiry = expiry
        self.strike = strike
        self.fwd = fwd
        self.riemann = riemann
        self.debug = debug
    
    def computeCallBlackPrice(self):
        self.q = (numpy.power(self.strike, 1.0 - self.beta))/(1.0 - self.beta)
        self.q0 = (numpy.power(self.fwd, 1.0 - self.beta))/(1.0 - self.beta)
        self.eta = numpy.abs(1.0/(2.0*(self.beta-1.0)))
        self.s_plus = numpy.arcsinh((self.nu*numpy.abs(self.q+self.q0))/self.alpha)
        self.s_plus_sinh = numpy.sinh(self.s_plus)
        self.s_minus = numpy.arcsinh((self.nu*numpy.abs(self.q-self.q0))/self.alpha)
        self.s_minus_sinh = numpy.sinh(self.s_minus)
        intrinsic_value = numpy.maximum(self.fwd-self.strike,0)
        integral = self.computeIntegral()
        price_pre_factor = (2.0/numpy.pi)*numpy.sqrt(self.fwd*self.strike)
        return intrinsic_value + price_pre_factor*integral 

    def computeIntegral(self):
        def kappa(s):
            sinh_s = numpy.sinh(s)
            sinh_s = sinh_s*sinh_s
            sq_rt = numpy.sqrt((sinh_s-self.s_minus_sinh*self.s_minus_sinh)/(self.s_plus_sinh*self.s_plus_sinh-sinh_s))
            ret = 2.0*numpy.arctan(sq_rt)
            return ret

        def psi(s):
            sinh_s = numpy.sinh(s)
            sinh_s = sinh_s*sinh_s
            sq_rt = numpy.sqrt((sinh_s-self.s_plus_sinh*self.s_plus_sinh)/(sinh_s-self.s_minus_sinh*self.s_minus_sinh))
            ret = 2.0*numpy.arctanh(sq_rt)
            return ret

        def KernelAnalytical(s):
            t = self.expiry*self.nu*self.nu
            s_sqr = s*s
            kernel_multiplic = numpy.sqrt(numpy.sinh(s)/s)*numpy.exp(-(s_sqr/(2.0*t))-(t/8.0))
            g = s*(numpy.cosh(s)/numpy.sinh(s))-1.0
            r = 1.0+(3/8)*((t*g)/s_sqr)-5.0*t*t*(-8.0*s_sqr+3.0*g*g+24.0*g)/(128.0*s_sqr*s_sqr)+35.0*t*t*t*(-40.0*s_sqr+3.0*g*g*g+24.0*g*g+120.0*g)/(1024.0*numpy.power(s,6.0))
            delta_r = numpy.exp(t/8.0)-(3072.0+384.0*t+24.0*t*t+t*t*t)/3072.0
            ret = kernel_multiplic*(r+delta_r)
            return ret

        def integrand_1(s):
            ret = (numpy.sin(self.eta*kappa(s))/numpy.sinh(s))*KernelAnalytical(s)
            return ret

        def integrand_2(s):
            ret = (numpy.exp(-self.eta*psi(s))/numpy.sinh(s))*KernelAnalytical(s)
            return ret
        
        try:
            integral1 = 0.0
            integral2 = 0.0
            '''
            @todo: when using Riemann, the ATM strike vol returned is 0, this 
            should be investigated. The same doesn't happen when using Numpy 
            integration functions, although the ATM result seems incorrect
            for some parameter combinations
            '''
#           Riemann integration approach
            if self.riemann:
                '''
                @todo: the following steps seem to be in line with the quadrature
                Numpy results, although it might be possible to further refine
                them to speed up the computation
                '''
                integral_1_step = 0.001        
                integral_2_step = 0.001
                '''
                @todo: the upper bound needs to be reviewed, in particular
                a method to find the right upper bound has to be developed
                as we might find that a lower bound would also do,
                increasing in this way the efficiency 
                '''
                integral_2_upper_bound = 25.0
                
                integration_index = self.s_minus
                
                while integration_index <= (self.s_plus - integral_1_step):
                    integral1 += (integrand_1(integration_index + integral_1_step/2.0))*integral_1_step
                    integration_index += integral_1_step

                integration_index = self.s_plus
                    
                if not self.beta == 0.5:
                    while integration_index <= integral_2_upper_bound:
                        integral2 += (integrand_2(integration_index + integral_2_step/2.0))*integral_2_step
                        integration_index += integral_2_step                  
            else:
                integral1 = integrate.fixed_quad(func=integrand_1,a=self.s_minus,b=self.s_plus)[0]
                '''
                Using quadrature for the following integral doesn't seem to work
                integral2 = integrate.quad(func=integrand_2,a=self.s_plus,b=inf)[0]
                '''                
                sub_integral2 = 1.0
                end_interval = self.s_plus
                while abs(sub_integral2)>abs(integral2)/1000.0:
                    init_interval = end_interval
                    '''
                    @todo: the following 0.001 should be reviewed,
                    possibly made coarser to speed up computation
                    '''
                    end_interval = init_interval+0.001*self.s_plus
                    sub_integral2 = integrate.fixed_quad(func=integrand_2, a=init_interval, b=end_interval, n=10)[0]
                    integral2 += sub_integral2
                    
            int2_pre_factor = numpy.sin(self.eta*numpy.pi)
            
            integral = integral1 + int2_pre_factor*integral2
                
            return integral
        except Exception as e:
            print 'Error in computeIntegral main routine: %s' %(e)
            