import numpy as np

class Clutch():
    def __init__(self, n_discs, P_in=5475, n_in=750, safety_factor=2, IR_fric=0.190, OR_fric=0.250, OR_spline=0.260, IR_spline=0.175, mu_fric=0.14, mu_spline=0.1, OR_hyd=0.270, IR_hyd=0.150, Oring_ef=0.9, clutch_length=0.165) -> None:
        '''
        n_discs: number of discs
        P_in: input power [kW]
        n_in: input speed [rpm]
        safety_factor: safety factor
        IR_fric: inner radius of friction material [m]
        OR_fric: outer radius of friction material [m]
        OR_spline: outer radius of spline [m]
        IR_spline: inner radius of spline [m]
        mu_fric: friction coefficient of friction material
        mu_spline: friction coefficient of spline
        OR_hyd: outer radius of hydraulic piston [m]
        IR_hyd: inner radius of hydraulic piston [m]
        Oring_ef: efficiency of Oring
        Clutch_length: length of clutch housing [m]
        '''
        self.n_discs = n_discs
        self.P_in = P_in
        self.n_in = n_in
        self.safety_factor = safety_factor
        self.IR_fric = IR_fric
        self.OR_fric = OR_fric
        self.OR_spline = OR_spline
        self.IR_spline = IR_spline
        self.mu_fric = mu_fric
        self.mu_spline = mu_spline
        self.OR_hyd = OR_hyd
        self.IR_hyd = IR_hyd
        self.Oring_ef = Oring_ef
        self.clutch_length = clutch_length

    def torque(self):
        rads = self.n_in * 2 * np.pi / 60
        self.T = self.P_in*1000 / rads
    
    def hydraulic_area(self):
        self.A = np.pi * (self.OR_hyd**2 - self.IR_hyd**2)
    
    def required_pressure(self, model=1):
        '''
        Computes the required normal force on the friction material
        model: 1 for uniform pressure model, 2 for uniform wear model
        Assumes first plate is mounted outside
        '''
        # Computes the torque
        self.torque()
        required_torque = self.T * self.safety_factor
        # Computes the hydraulic area
        self.hydraulic_area()


        if model == 1: # Uniform pressure model
            r_f = 2/3 * (self.IR_fric**3 - self.OR_fric**3) / (self.IR_fric**2 - self.OR_fric**2) 
        if model == 2: # Uniform wear model
            r_f = (self.IR_fric + self.OR_fric) / 2

        inner_coef = (1 - r_f / self.IR_spline * self.mu_fric * self.mu_spline) / (1 + r_f / self.IR_spline * self.mu_fric * self.mu_spline)
        outer_coef = (1 - r_f / self.OR_spline * self.mu_fric * self.mu_spline) / (1 + r_f / self.OR_spline * self.mu_fric * self.mu_spline)

        torque_conv = r_f * self.mu_fric

        denom1 = self.Oring_ef * outer_coef * torque_conv
        
        denom2 = 0
        a=0
        b=0
        for i in range(self.n_discs):
            denom2_i = 0
            if i == 0:
                denom2_i = 1
            
            elif (i % 2) == 1: # For uneven numbers
                a += 1
                denom2_i = inner_coef ** (a) * outer_coef ** (a-1)
            
            elif (i % 2) == 0: # For even numbers
                b += 1
                denom2_i = inner_coef ** (b) * outer_coef ** (b)
            denom2 = denom2 + denom2_i

        F0 = required_torque / (denom1 * denom2)

        P = F0 / self.A

        P = P / 10**5 # Convert to Bar

        return P
    
    def disc_width(self):
        disc_width = self.clutch_length/self.n_discs
        disc_width = disc_width * 1000 # Convert to mm
        return disc_width




def find_gear_pair(i_goal = 292/750, max_teeth=150, min_teeth=20, min_beta=5, max_beta=20, beta_steps=5, distance=0.725, max_profile_shift=2):
    # All combinations of teeth
    z1 = np.arange(min_teeth,max_teeth+1,1)
    z2 = np.rint(z1 / i_goal)

    # Trim z2 vector such that the highest number is under 150
    z2 = z2[z2 <= max_teeth]
    # Trim z1 vector to match the length of z2
    z1 = z1[:len(z2)]

    alpha_w = distance
    possible_beta = np.arange(min_beta, max_beta+beta_steps, beta_steps)*np.pi/180
    possible_module = np.array([6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18])*0.001
    alpha_n = 20 * np.pi / 180

    possible_gear_pairs = np.empty((0,5))

    for j, beta in enumerate(possible_beta):
        profile_shift = np.zeros((len(possible_module), len(z1)))
        for i, m_t in enumerate(possible_module):
            for n in range(len(z1)):
                a = m_t * (z1[n] + z2[n]) / 2
                alpha_t = np.arctan(np.tan(alpha_n) / np.cos(beta))
                invalpha_t = np.tan(alpha_t) - alpha_t
                alpha_wt = np.arccos(a * np.cos(alpha_t) / alpha_w)
                invalpha_wt = np.tan(alpha_wt) - alpha_wt
                profile_shift[i, n] = (z1[n] + z2[n]) / 2 * (invalpha_wt - invalpha_t) / np.tan(alpha_t)
        profile_shift[profile_shift < 0] = np.nan
        profile_shift[profile_shift > max_profile_shift] = np.nan

        for n in range(len(profile_shift[0,:])):
            for i in range(len(profile_shift[:,0])):
                if not np.isnan(profile_shift[i,n]):
                    possible = [z1[n], z2[n], possible_module[i], profile_shift[i,n], beta]
                    possible_gear_pairs = np.append(possible_gear_pairs, [possible], axis=0)

    zero_columns = np.where(possible_gear_pairs[0] == 0)[0]
    possible_gear_pairs = np.delete(possible_gear_pairs, zero_columns, axis=1)

    return possible_gear_pairs
    

class GearBox():
    def __init__(self, z_1, z_2, m_n, profile_shift, beta, P_in=5475, n_in = 750, E=210000 , b=0.375, a_w=0.725):
        '''
        n_1: number of teeth on input gear
        n_2: number of teeth on output gear
        m_n: module of gears
        profile_shift: profile shift of gears
        beta: helix angle
        P_in: input power [kW]
        n_in: input speed [rpm]
        E: Young's modulus [MPa]
        b: face width [m]
        a_w: axel distance [m]
        '''
        self.z_1 = z_1
        self.z_2 = z_2
        self.m_n = m_n
        self.alpha_n = 20 * np.pi / 180
        self.profile_shift = profile_shift
        if self.profile_shift > 1:
            self.x1 = 0.5
            self.x2 = self.profile_shift - 0.5
        else:
            self.x1 = self.profile_shift * 2/3
            self.x2 = self.profile_shift * 1/3
        self.x1 = 0.5
        self.x2 = self.profile_shift - 0.5
        self.beta = beta
        self.b = b
        self.a_w = a_w
        self.P_in = P_in
        self.n_in = n_in
        self.E = E

        self.run_all()
        
    def transverse_module(self):
        self.m_t = self.m_n * np.cos(self.beta)
        return self.m_t
    
    def transverse_pressure_angle(self):
        self.alpha_t = np.arctan(np.tan(self.alpha_n) / np.cos(self.beta))
        return self.alpha_t
        
    def reference_circle_radius(self):
        self.r_1 = self.m_t * self.z_1 / 2
        self.r_2 = self.m_t * self.z_2 / 2
        return self.r_1, self.r_2

    def reference_circle_pitch(self):
        self.p_t = np.pi * self.m_t
        return self.p_t

    def base_circle_radius(self):
        self.r_b1 = self.r_1 * np.cos(self.alpha_t)
        self.r_b2 = self.r_2 * np.cos(self.alpha_t)
        return self.r_b1, self.r_b2

    def base_circle_pitch(self):
        self.p_et = self.p_t * np.cos(self.alpha_t)
        return self.p_et

    def base_helix_angle(self):
        self.beta_b = np.arcsin(np.sin(self.beta) * np.cos(self.alpha_n))
        return self.beta_b
    
    def unmodified_tip_circle_radius(self):
        self.r_a1 = self.r_1 + (1 + self.x1) * self.m_n
        self.r_a2 = self.r_2 + (1 + self.x2) * self.m_n
        return self.r_a1, self.r_a2
    
    def root_circle_radius(self):
        self.r_f1 = self.r_1 - (1.25 - self.x1) * self.m_n
        self.r_f2 = self.r_2 - (1.25 - self.x2) * self.m_n
        return self.r_f1, self.r_f2
    
    def reference_center_distance(self):
        self.a = self.m_t * (self.z_1 + self.z_2) / 2
        return self.a
    
    def working_transverse_pressure_angle(self):
        self.alpha_wt = np.arccos(self.a * np.cos(self.alpha_t) / self.a_w)
        return self.alpha_wt
    
    def recommended_tip_circle_radius(self):
        t = (self.a - self.a_w) / self.m_n + self.x1 + self.x2-0.1
        if t < 0:
            t = 0

        self.r_a1 = self.r_1 + (1 + self.x1 - t) * self.m_n
        self.r_a2 = self.r_2 + (1 + self.x2 - t) * self.m_n
        return self.r_a1, self.r_a2
    
    def transverse_contact_ratio(self):
        self.epsilon_alpha = (np.sqrt(self.r_a1**2 - self.r_b1**2) + np.sqrt(self.r_a2**2 - self.r_b2**2) - (self.r_b1 + self.r_b2) * np.tan(self.alpha_wt)) / (np.pi * self.m_t * np.cos(self.alpha_t))
        return self.epsilon_alpha
    
    def overlap_ratio(self):
        self.epsilon_beta = self.b * np.tan(self.beta_b) / self.p_et
        return self.epsilon_beta
    
    def total_contact_ratio(self):
        self.epsilon_gamma = self.epsilon_alpha + self.epsilon_beta
        return self.epsilon_gamma

    def shaft_torque(self):
        rads = self.n_in * 2 * np.pi / 60
        self.T_1 = self.P_in*1000 / rads
        return self.T_1
    
    def tramsverse_gear_force(self):
        self.F_bt = self.T_1 / self.r_b1
        return self.F_bt
    
    def gear_normal_force(self):
        self.F_bn = self.F_bt / np.cos(self.beta_b)
        return self.F_bn
    
    def axial_gear_force(self):
        self.F_a = self.F_bn * np.sin(self.beta_b)
        return self.F_a
    
    def transverse_tangential_gear_force(self):
        self.F_t = self.F_bt * np.cos(self.alpha_wt)
        return self.F_t
    
    def transverse_radial_gear_force(self):
        self.F_r = self.F_bt * np.sin(self.alpha_wt)
        return self.F_r

    def run_all(self):
        self.transverse_module()
        self.transverse_pressure_angle()
        self.reference_circle_radius()
        self.reference_circle_pitch()
        self.base_circle_radius()
        self.base_circle_pitch()
        self.base_helix_angle()
        self.unmodified_tip_circle_radius()
        self.root_circle_radius()
        self.reference_center_distance()
        self.working_transverse_pressure_angle()
        self.recommended_tip_circle_radius()
        self.transverse_contact_ratio()
        self.overlap_ratio()
        self.total_contact_ratio()
        self.shaft_torque()
        self.tramsverse_gear_force()
        self.gear_normal_force()
        self.axial_gear_force()
        self.transverse_tangential_gear_force()
        self.transverse_radial_gear_force()

    def b_correction(self, update_b=False):
        '''
        Finds the correction of b such that the total contact ratio is an int.
        '''
        self.nearest_int = round(self.epsilon_gamma)
        self.distance_to_int = (self.nearest_int - self.epsilon_gamma)
        correction = self.distance_to_int * self.p_et / np.tan(self.beta_b)
        self.b_corrected = self.b + correction
        if update_b:
            self.b = self.b_corrected
            self.run_all()
        return self.b_corrected
    
    def saftey_factor(self):

    
        #Application factor, KA
        K_A = 1 # Side 242

        # Dynamic factor, KV - Quality class 6
        f_F = 1 # Side 243
        K = 32 # Side 243
        v = 2*np.pi * self.r_1 * self.n_in / 60
        u = self.z_2 / self.z_1 
        K_V = 1 + f_F * K * self.z_1 * v * np.sqrt(u**2/(1 + u**2))*10e-6
        w_t = self.F_bt / self.b * K_A * K_V # Side 243
        
        # Axial load distribution factor, K_Hbeta, K_Fbeta
        K_beta = 1.15 # Side 245 width skulle gerne være 315 eller derover
        if w_t < 100: # Side 245
            f_w = 1.6
        elif 100 <= w_t < 200:
            f_w = 1.45
        elif 200 <= w_t < 250:
            f_w = 1.3
        elif 250 <= w_t < 315:
            f_w = 1.15
        else:
            f_w = 1.0
        
        f_p = 1 # Side 245 - Steel
        K_Fbeta = 1 + (K_beta - 1) * f_w * f_p

        K_Hbeta = K_Fbeta**1.39 # Side 244

        # Transverse load distribution factor, K_Halpha, K_Falpha
        f_pe = 15
        y_p = 2.5*10**(-6)
        w_t = self.F_bt / self.b * K_A * K_V # Side 243

        if self.epsilon_gamma <= 2:
            K_Halpha = self.epsilon_gamma / 2 * (0.9 + 0.4 * 20 * (f_pe - y_p) / (w_t * K_Fbeta) ) # Side 246
        else:
            K_Halpha = 0.9+0.4* np.sqrt(2*(self.epsilon_gamma-1)/self.epsilon_gamma) * 20 * (f_pe - y_p) / (w_t * K_Fbeta) # Side 247 

        if K_Halpha < 1:
            K_Halpha = 1
        

        # Pitting
        #Zone Factor
        Z_H = np.sqrt((2 * np.cos(self.beta_b))/ np.cos(self.alpha_t)**2 * self.alpha_t * np.tan(self.alpha_wt)) # Side 249

        #Elasitcity factor
        Z_E = np.sqrt(0.175*self.E) # side 249

        # Contact ratio factor
        if self.epsilon_beta < 1:
            Z_epsilon = np.sqrt((4-self.epsilon_alpha)/ 3 * (1 - self.epsilon_beta) + self.epsilon_beta / self.epsilon_alpha) # Side 249
        else:
            Z_epsilon = np.sqrt(1 / self.epsilon_alpha) # Side 249

        #Helix angle factor
        Z_beta = np.sqrt(np.cos(self.beta))

        # Life factor
        N_L = 10**6 # Side 250
        Z_N = (10**9 / N_L)**0.0057 # Side 250 - Hardenend steel with some pitting allowed

        # Lubrication factor
        sigma_Hlim = 1000 # Side 251 - Gæt
        C_ZL = (sigma_Hlim - 850) / 350 * 0.08 + 0.83

        v_40 = 100

        Z_L = C_ZL + (1 * (1-C_ZL)) / (1.2 + 134 / v_40)

        # Roughness factor
        R_z1 = 3 # ISO Standard
        R_z2 = 3
        R_z100 = (R_z1 + R_z2) / 2 * np.cbrt(100/self.a)
        C_ZR = 0.12 + (1000 - sigma_Hlim) / 5000

        Z_R = (3/R_z100)**C_ZR

        # Speed factor
        C_ZV = 0.85 + (sigma_Hlim - 850) / 350 * 0.08
        Z_V = C_ZV + (2*(1 - C_ZV))/np.sqrt(0.8 + 32/v)

        # Work hardening factor
        Z_W = 1 # ISO standard

        Z_X = 1 # side 248

        sigma_H01 = Z_H * Z_E * Z_epsilon * Z_beta * np.sqrt(self.F_bt/(self.r_1*2*self.b) * (u + 1)/u) # Side 248
        sigma_H02 = Z_H * Z_E * Z_epsilon * Z_beta * np.sqrt(self.F_bt/(self.r_2*2*self.b) * (u + 1)/u)

        S_H1 = sigma_Hlim * Z_N / sigma_H01 * (Z_L * Z_R * Z_V * Z_W * Z_X) / np.sqrt(K_A * K_V * K_Hbeta * K_Halpha)
        S_H2 = sigma_Hlim * Z_N / sigma_H02 * (Z_L * Z_R * Z_V * Z_W * Z_X) / np.sqrt(K_A * K_V * K_Hbeta * K_Halpha)

        # Stress
        # Tooth form factor
        x1 = round(abs(self.x1),1)
        x2 = round(abs(self.x2),1)
        if self.z_1 < 55:
            if x1 == 0:
                Y_FS1 = 4.66
            elif x1 == 0.1:
                Y_FS1 = 4.55
            elif x1 == 0.2:
                Y_FS1 = 4.45
            elif x1 == 0.3:
                Y_FS1 = 4.38
            elif x1 == 0.4:
                Y_FS1 = 4.33
            elif x1 == 0.5:
                Y_FS1 = 4.3
            elif x1 == 0.6:
                Y_FS1 = 4.27
            elif x1 == 0.7:
                Y_FS1 = 4.25
            elif x1 == 0.8:
                Y_FS1 = 4.22
            elif x1 == 0.9:
                Y_FS1 = 4.18
            else:
                Y_FS1 = 4.12
        if self.z_1 >= 55:
            if x1 == 0:
                Y_FS1 = 4.26
            elif x1 == 0.1:
                Y_FS1 = 4.28
            elif x1 == 0.2:
                Y_FS1 = 4.3
            elif x1 == 0.3:
                Y_FS1 = 4.33
            elif x1 == 0.4:
                Y_FS1 = 4.35
            elif x1 == 0.5:
                Y_FS1 = 4.38
            elif x1 == 0.6:
                Y_FS1 = 4.41
            elif x1 == 0.7:
                Y_FS1 = 4.43
            elif x1 == 0.8:
                Y_FS1 = 4.43
            elif x1 == 0.9:
                Y_FS1 = 4.43
            else:
                Y_FS1 = 4.41  
        
        if self.z_2 < 55:
            if x2 == 0:
                Y_FS2 = 4.66
            elif x2 == 0.1:
                Y_FS2 = 4.55
            elif x2 == 0.2:
                Y_FS2 = 4.45
            elif x2 == 0.3:
                Y_FS2 = 4.38
            elif x2 == 0.4:
                Y_FS2 = 4.33
            elif x2 == 0.5:
                Y_FS2 = 4.3
            elif x2 == 0.6:
                Y_FS2 = 4.27
            elif x2 == 0.7:
                Y_FS2 = 4.25
            elif x2 == 0.8:
                Y_FS2 = 4.22
            elif x2 == 0.9:
                Y_FS2 = 4.18
            else:
                Y_FS2 = 4.12
        if self.z_2 >= 55:
            if x2 == 0:
                Y_FS2 = 4.26
            elif x2 == 0.1:
                Y_FS2 = 4.28
            elif x2 == 0.2:
                Y_FS2 = 4.3
            elif x2 == 0.3:
                Y_FS2 = 4.33
            elif x2 == 0.4:
                Y_FS2 = 4.35
            elif x2 == 0.5:
                Y_FS2 = 4.38
            elif x2 == 0.6:
                Y_FS2 = 4.41
            elif x2 == 0.7:
                Y_FS2 = 4.43
            elif x2 == 0.8:
                Y_FS2 = 4.43
            elif x2 == 0.9:
                Y_FS2 = 4.43
            else:
                Y_FS2 = 4.41

        # Helix angle factor
        Y_beta = 1 - self.epsilon_beta * self.beta/(2*np.pi/3)
        y_betamin = 1 - 0.25*self.beta
        if y_betamin < 0.75:
            y_betamin = 0.75
        if Y_beta < y_betamin:
            Y_beta = y_betamin

        # Life factor
        N_L = 10**4
        Y_NT = (3*10**6 / N_L)**0.12 # 254

        # Relative notch sensitivity
        Y_beta = 1

        #Relative surface condition factor
        Y_R = 1.674 - 0.529*(R_z1+1)**0.1 # 255 - 12.71

        # Size factor
        if self.m_n < 0.003:
            Y_X = 1.05 - 0.01 * self.m_n
        else:
            Y_X = 0.75 # side 255

        Y_epsilon = 0.25 + 0.75 / self.epsilon_alpha # side 247

        sigma_F01 = self.F_bt / (self.b* self.m_n*1000) * Y_FS1 * Y_epsilon * Y_beta # side 252
        sigma_F02 = self.F_bt / (self.b* self.m_n*1000) * Y_FS2 * Y_epsilon * Y_beta

        sigma_Flim = 300 #ISO - Figure 1

        S_F1 = (sigma_Flim * Y_NT) / sigma_F01 * (Y_beta * Y_R * Y_X) / (K_A * K_V * K_Fbeta * K_Halpha)
        S_F2 = (sigma_Flim * Y_NT) / sigma_F02 * (Y_beta * Y_R * Y_X) / (K_A * K_V * K_Fbeta * K_Halpha)

        return S_H1, S_F1, S_F2



if __name__ == "__main__":
    if False: # Question 2
        discs = np.arange(4, 50, 1)
        # discs = [6]
        for i, disc in enumerate(discs):
            clutch = Clutch(n_discs=disc)
            # print(disc, clutch.required_pressure(model=2))
            print(disc, clutch.disc_width())

    if False: # Question 3
        possible_gears  = find_gear_pair()
        print(possible_gears)

    if False: # Question 4
        possible_gears  = find_gear_pair()

        gears = []

        for i, gear in enumerate(possible_gears):
            gears.append(GearBox(gear[0], gear[1], gear[2], gear[3], gear[4]))

        for i, gear in enumerate(gears):
            if gear.beta == 5 * np.pi / 180:
                print(gear.z_1, gear.m_n)
                # print(gear.z_1, gear.profile_shift)     

    if False: # question 5
        possible_gears  = find_gear_pair()        

        gears = []

        for i, gear in enumerate(possible_gears):
            gears.append(GearBox(gear[0], gear[1], gear[2], gear[3], gear[4]))

        # Remove gears with a negative epsilon_alpha
        y = [s for s in gears if s.epsilon_alpha > 0]

        for i, gear in enumerate(y):
            if gear.beta == 15 * np.pi / 180:
                print(gear.z_1, gear.epsilon_alpha)
                # print(gear.z_1, gear.epsilon_gamma)
                # print(gear.z_1, gear.b_correction())

    if False: # Question 6
        possible_gears  = find_gear_pair()        

        gears = []

        for i, gear in enumerate(possible_gears):
            gears.append(GearBox(gear[0], gear[1], gear[2], gear[3], gear[4]))

        # Remove gears with a negative epsilon_alpha
        y = [s for s in gears if s.epsilon_alpha > 0]

        for i, gear in enumerate(y):
            if gear.beta == 15 * np.pi / 180:
                print(gear.z_1, gear.saftey_factor())

    if False: # Geat selection
        possible_gears  = find_gear_pair()        

        gears = []

        for i, gear in enumerate(possible_gears):
            gears.append(GearBox(gear[0], gear[1], gear[2], gear[3], gear[4]))

        # Remove gears with a negative epsilon_alpha
        y = [s for s in gears if s.epsilon_alpha > 0]

        for i, gear in enumerate(y):
            if gear.beta == 10 * np.pi / 180:
                if gear.z_1 == 27:
                    print(gear.z_1, gear.z_2, gear.m_n, gear.profile_shift, gear.epsilon_gamma, gear.saftey_factor())

    if True: # Forces
        possible_gears  = find_gear_pair()        

        gears = []

        for i, gear in enumerate(possible_gears):
            gears.append(GearBox(gear[0], gear[1], gear[2], gear[3], gear[4]))

        # Remove gears with a negative epsilon_alpha
        y = [s for s in gears if s.epsilon_alpha > 0]

        for i, gear in enumerate(y):
            if gear.beta == 10 * np.pi / 180:
                if gear.z_1 == 27:
                    print(gear.F_t, gear.F_r, gear.F_a, gear.r_1)

    if False: # Extra
        possible_gears  = find_gear_pair()        

        gears = []

        for i, gear in enumerate(possible_gears):
            gears.append(GearBox(gear[0], gear[1], gear[2], gear[3], gear[4]))

        # Remove gears with a negative epsilon_alpha
        y = [s for s in gears if s.epsilon_alpha > 0]

        for i, gear in enumerate(y):
            # print(gear.z_1, gear.x1)
            print(gear.z_1, gear.x2)
