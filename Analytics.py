import random
import math
import numpy
import scipy
from scipy . stats import norm
numpy . random . seed (0)
# Seed the random number generator
random . seed(0)

class Analytics():
    
    def __init__(self):
        self._debug = False
        self.small_number = 1e-06
        
    def letSmallNumber(self):
        return self.small_number
    
    def ComputeFirstDerivative(self, f_plus, f_minus):
        try:
            derivative = (f_plus - f_minus) / (2.0*self.small_number)
            return derivative
        except Exception as e:
            print 'Error in ComputeFirstDerivative: %s' %(e)

    def ComputeSecondDerivative(self, f, f_plus, f_minus):
        try:
            derivative = (f_plus - 2.0*f + f_minus) / (self.small_number*self.small_number)
            return derivative
        except Exception as e:
            print 'Error in ComputeSecondDerivative: %s' %(e)

    def black (self, F_0 , y, expiry , vol , isCall ):
        option_value = 0
        if expiry * vol == 0.0:
            if isCall :
                option_value = max(F_0 - y, 0.0)
            else :
                option_value = max(y - F_0 , 0.0)
        else :
            d1 = self.dPlusBlack (F_0 = F_0 , y = y, expiry = expiry , vol = vol)
            d2 = self.dMinusBlack (F_0 = F_0 , y = y, expiry = expiry , vol = vol)
            if isCall :
                option_value = (F_0 * norm .cdf(d1) - y *norm .cdf (d2 ))
            else :
                option_value = (y * norm .cdf(-d2) - F_0 * norm .cdf (-d1 ))
        return option_value
    
    def dPlusBlack (self, F_0 , y, expiry , vol ):
        d_plus = (( math .log (F_0 / y) + 0.5 * vol * vol * expiry ) / vol / math . sqrt ( expiry ))
        return d_plus
    
    def dMinusBlack (self, F_0 , y, expiry , vol ):
        d_minus = (self.dPlusBlack (F_0 = F_0 , y = y, expiry = expiry , vol = vol) - vol * math . sqrt ( expiry ))
        return d_minus
            
    def ComputeStandardSABRBlackVol(self, fwd, strike, maturity, alpha, beta, rho, nu):
        try:
            one_minus_beta = scipy.power(1.0-beta, 2.0)
            if fwd == strike:
                term_1 = alpha/scipy.power(fwd, (1.0-beta))           
                term_2 = (one_minus_beta/24.0)*alpha*alpha/(scipy.power(fwd, (2.0-2.0*beta)))
                term_3 = 0.25*(rho*beta*nu*alpha/scipy.power(fwd, (1.0-beta)))
                term_4 = ((2.0-3.0*rho*rho)/24.0)*nu*nu
                sigma = term_1*(1.0+(term_2+term_3+term_4)*maturity)
            else:
                f_k_pow = scipy.power(fwd*strike, (1.0-beta)/2.0)
                log_f_k = scipy.log(fwd/strike)
                z = (nu/alpha)*f_k_pow*log_f_k
                x = scipy.log((scipy.sqrt(1.0-2.0*rho*z+z*z)+z-rho)/(1.0-rho))
                term_1 = alpha/(f_k_pow*(1.0+(one_minus_beta/24.0)*log_f_k*log_f_k+(one_minus_beta*one_minus_beta/1920.0)*log_f_k*log_f_k*log_f_k*log_f_k))
                term_2 = z/x
                term_3 = (one_minus_beta/24.0)*alpha*alpha/(f_k_pow*f_k_pow)
                term_4 = 0.25*(rho*beta*nu*alpha/f_k_pow)
                term_5 = ((2.0-3.0*rho*rho)/24.0)*nu*nu
                term_6 = 1.0+(term_3+term_4+term_5)*maturity
                sigma = term_1*term_2*term_6
            if self._debug:
                print 'SABR Black Vol: %s' %(sigma)
            return sigma
        except Exception as e:
            print 'Error in ComputeStandardSABRBlackVol: %s' %(e)
            
    def drawTwoRandomNumbers(self, rho):
        '''
        Draw a pair of correlated random numbers .
        @var rho : SABR Rho
        '''
        try:
            rand_list = []
            z1 = numpy.random.normal()
            y1 = numpy.random.normal()
            rand_list.append (z1)
            term1 = z1 * rho
            term2 = (y1 *
                     math . pow ((1.0 - math . pow(rho , 2.0)) , 0.5))
            x2 = term1 + term2
            rand_list.append (x2)
            return rand_list
        except Exception as e:
            print 'Error in drawTwoRandomNumbers: %s' %(e)
   
    def computeSABRMonteCarloEuler(self, no_of_sim , no_of_steps ,
    maturity , F_0 , alpha_0 , beta , rho , nu ):
        '''
        Monte Carlo SABR using Euler scheme .
        @var no_of_sim : Monte Carlo paths
        @var no_of_steps : discretization steps required
        to reach the option maturity date
        @var maturity : option maturity (in years )
        @var F_0 : forward interest rate
        @var alpha_0 : SABR Alpha at t=0
        @var beta : SABR Beta
        @var rho : SABR Rho
        @var nu : SABR Nu
        '''
        try:
            # Step length in years
            dt = float(maturity) / float(no_of_steps)
            dt_sqrt = math.sqrt(dt)
            simulated_forwards = []
            no_of_sim_counter = 0
            while no_of_sim_counter < no_of_sim :
                F_t = F_0
                alpha_t = alpha_0
                no_of_steps_counter = 1
                while no_of_steps_counter <= no_of_steps :                
                    # Zero absorbing boundary used for all the beta
                    # choices except beta = 0 and beta = 1

                    if ((beta > 0 and beta < 1) and F_t <= 0):
                        F_t = 0
                        no_of_steps_counter = no_of_steps + 1
                    else:
                        #Generate two correlated
                        # random numbers
                        rand = self.drawTwoRandomNumbers(rho)
                        w_F = dt_sqrt * rand[0]
                        F_b = math.pow(abs(F_t), beta)
                        F_t = F_t + alpha_t * F_b * w_F
                        w_a = dt_sqrt * rand[1]
    #                   alpha_t = (alpha_t + nu * alpha_t * w_a)
                        # Use the actual solution for the vol SDE
                        alpha_t = alpha_t * math.exp(nu * w_a - 0.5 * nu * nu * dt)
                        
                    no_of_steps_counter += 1
                # At the end of each path , we store the forward
                # interest rate in a list
                simulated_forwards.append(F_t)
                no_of_sim_counter = no_of_sim_counter + 1
                
            return simulated_forwards
        except Exception as e:
            print 'Error in computeSABRMonteCarloEuler: %s' %(e)
                        
if __name__ == '__main__':
    strikes = [0.01, 0.02, 0.03, 0.04, 0.05]
    maturity = 19.5
    fwd = 0.021820668

    beta = 0.10162
    rho = 0.37817
    alpha = 0.00982
    nu = 0.40425
    
    no_of_sim = 100000
    no_of_steps = 4
     
    analytics = Analytics()

    simulated_forwards = analytics.computeSABRMonteCarloEuler(no_of_sim, no_of_steps, maturity, fwd, alpha, beta, rho, nu)
    print simulated_forwards