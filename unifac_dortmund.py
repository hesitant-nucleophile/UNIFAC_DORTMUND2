import pandas as pd
import numpy as np
from rdkit import Chem
from ugropy import Groups
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

# load from excel
file_path1 = os.path.join(script_dir, "interaction_parameters.xlsx")
df_ip = pd.read_excel(file_path1)
df_ip = df_ip.fillna(0) # replaces nan with 0

file_path2 = os.path.join(script_dir, "groupe_surfaces_and_volumes.xlsx")
df_sv=pd.read_excel(file_path2)


class Unifac_Dortmund:
    def __init__(self,smiles_lst,mol_lst,T):
        self.smiles_lst=smiles_lst
        self.mol_lst=mol_lst
        self.T=T
        self.list_groups, self.unique_group_no = self.format_data()

        # error messages
        if len(mol_lst) != len(smiles_lst):
            raise ValueError(f"phi_lst (len={len(mol_lst)}) and smiles_lst (len={len(smiles_lst)}) must have same length.")
        
        if sum(mol_lst)!=1:
            raise ValueError(f"molefractions must sum to unity, but where of sum {sum(mol_lst)}")
    
    # smiles to groupes
    def disect_smiles(self,smiles):
        molecule = Groups(smiles, identifier_type="smiles")
        unifac_counts = molecule.dortmund.subgroups
        return unifac_counts
    
    # making list of dataframes with group number, group type and number of that groupe
    def format_data(self):
        smiles_lst=self.smiles_lst
        data = {"Group_No": pd.Series(dtype=int),
        "Group": pd.Series(dtype=str),
        "number": pd.Series(dtype=int)}

        list_groups=[]

        for i in range(len(smiles_lst)):
            df = pd.DataFrame(data)
            smiles=smiles_lst[i]
            group_library=self.disect_smiles(smiles)

            for j in range(len(group_library)):
                group = list(group_library.keys())[j]
                index = df_sv.index[df_sv["Subgroup Name"] == group][0]
                group_no=df_sv.loc[index,"No."] 
                number=group_library[group]

                df.loc[len(df),"Group_No"]=group_no
                df.loc[len(df)-1,"Group"]=group
                df.loc[len(df)-1,"number"]=number
            list_groups.append(df)

        # unique groups
        df2=list_groups[0]
        df2=df2["Group_No"]
        for i in range(1,len(list_groups)):
            df3=list_groups[i]["Group_No"]
            df2 = pd.concat([df2, df3], ignore_index=True)
        lst =df2.tolist() 

        unique_group_no = list(dict.fromkeys(lst))

        return list_groups,unique_group_no 
    
    # get main group number from subgroup number
    def get_main_group(self, group_no):
        idx = df_sv.index[df_sv["No."] == group_no][0]
        return df_sv.loc[idx, "Main Group No."]
    
    #============================fractions based on whole mixture===============================
    def cap_x(self,group_no):
        sum_j1=0
        sum_j2=0
        list_groups=self.list_groups
        
        mol_lst=self.mol_lst
        #sum over all chemical species in mixture
        for j in range(len(list_groups)):
            df=list_groups[j] #select species df
            
            try:
                index = df.index[df["Group_No"] == group_no][0] # obs, problem her hvis en forbindelse er i et molekyle og ikke det andet, lige nu er den bare defineret sum nul i dataen, man kan evt, sige if none, continue uden at definere den sim 0 
            except:
                continue
            #if index==none   continue 
            v_mj=df.loc[index,"number"]
            
            x_j=mol_lst[j]

            sum_j1+=v_mj*x_j
             
        for j in range(len(list_groups)): #sum over all chemical species in mixture
            df=list_groups[j]
            x_j=mol_lst[j]
            sum_n=0
            for n in range(len(df)): # sum over all groups in molecule
                
                v_nj=df.loc[n,"number"]
                sum_n+=v_nj*x_j
            sum_j2+=sum_n

        return sum_j1/sum_j2
    
    def cap_theta(self,group_no): 
        unique_group_no=self.unique_group_no

        index = df_sv.index[df_sv["No."] == group_no][0]
        Q_m=df_sv.loc[index,"Q"] #get from group contribution table
        
        cap_x_m=self.cap_x(group_no)
        
        sum_n=0
        for n in range(len(unique_group_no)): #sum over all groups present in mixture
            group_no2=unique_group_no[n]

            index = df_sv.index[df_sv["No."] == group_no2][0]
            Q_n=df_sv.loc[index,"Q"] #get from group contribution table
            
            cap_x_n= self.cap_x(group_no2)

            sum_n+=Q_n*cap_x_n
        cap_theta_m=cap_x_m*Q_m/sum_n
        return cap_theta_m

    #============================fractions based on molecule=============================
    def cap_x_i(self,group_no,species_no):
        list_groups=self.list_groups
        
        df=list_groups[species_no] #select species df
            
        try:
            index = df.index[df["Group_No"] == group_no][0] # obs, problem her hvis en forbindelse er i et molekyle og ikke det andet, lige nu er den bare defineret sum nul i dataen, man kan evt, sige if none, continue uden at definere den sim 0 
            v_mj=df.loc[index,"number"]
        except:
            v_mj=0
    
        sum_n=0
        for n in range(len(df)): # sum over all groups in molecule 
            v_nj=df.loc[n,"number"]
            sum_n+=v_nj

        return v_mj/sum_n
    

    def cap_theta_i(self,group_no, species_no):# pillede ved den her
        list_groups=self.list_groups

        index = df_sv.index[df_sv["No."] == group_no][0]
        Q_m=df_sv.loc[index,"Q"] #get from group contribution table
        
        df=list_groups[species_no]
        index2 = df.index[df["Group_No"] == group_no][0]
        #nu_m=df.loc[index2, "number"]
        cap_x_m=self.cap_x_i(group_no,species_no)
        
        sum_n=0
        for n in range(len(df)): #sum over all groups present in molecule
            group_no2=df.loc[n,"Group_No"]

            index = df_sv.index[df_sv["No."] == group_no2][0]
            Q_n=df_sv.loc[index,"Q"] #get from group contribution table
            
            index2 = df.index[df["Group_No"] == group_no2][0]
            #nu_n=df.loc[index2, "number"]
            cap_x_n= self.cap_x_i(group_no2,species_no)
            #sum_n+=nu_n*Q_n
            sum_n+=cap_x_n*Q_n
        #cap_theta_m_i=Q_m*nu_m/sum_n
        cap_theta_m_i=Q_m*cap_x_m/sum_n
        return cap_theta_m_i

    #=======================================ln gamma_k=========================================
    def ln_cap_gamma(self,group_no_1):
        T=self.T
        unique_group_no=self.unique_group_no

        index = df_sv.index[df_sv["No."] == group_no_1][0]
        cap_q_k=df_sv.loc[index,"Q"] #get from group contribution table
        
        sum_m_2=0
        sum_m_1=0
        #==========================sum m 1:===========================
        for m in range(len(unique_group_no)): # sum over all gropues in the mixture
            
            group_no_2=unique_group_no[m]
            
            cap_theta_m=self.cap_theta(group_no_2)
            # convert subgroup number to main group number for interaction parameters
            mg1 = self.get_main_group(group_no_1) #k
            mg2 = self.get_main_group(group_no_2) #m 

            
            if mg1==mg2:
                a_mk,a_km=0,0 # self interaction
            else:
                try:
                    # have to do try except as ij is in same row as ji and does not have thei own row
                    try:
                        index=df_ip.index[(df_ip["i"] == mg2) & (df_ip["j"] == mg1)][0]
                        Amk,Bmk,Cmk=df_ip.loc[index,"Aij"],df_ip.loc[index,"Bij"],df_ip.loc[index,"Cij"]
                        Akm,Bkm,Ckm=df_ip.loc[index,"Aji"],df_ip.loc[index,"Bji"],df_ip.loc[index,"Cji"]
                    except:
                        index=df_ip.index[(df_ip["i"] == mg1) & (df_ip["j"] == mg2)][0]
                        Amk,Bmk,Cmk=df_ip.loc[index,"Aij"],df_ip.loc[index,"Bij"],df_ip.loc[index,"Cij"]
                        Akm,Bkm,Ckm=df_ip.loc[index,"Aji"],df_ip.loc[index,"Bji"],df_ip.loc[index,"Cji"]
                except:
                    print(f"Interaction parameters not found for groups {mg1} and {mg2}")
                    return None
                a_mk=Amk+Bmk*T+Cmk*T**2
                a_km=Akm+Bkm*T+Ckm*T**2
            cap_psi_mk=np.exp(-a_mk/T)
            cap_psi_km=np.exp(-a_km/T)
            
            sum_m_1+=cap_theta_m*cap_psi_mk

            
            sum_n_1=0
            #==============================sum m 2:=========================
            #sum n:
            for n in range(len(unique_group_no)):
                group_no_3=unique_group_no[n]
                cap_theta_n=self.cap_theta(group_no_3)
                # convert subgroup number to main group number for interaction parameters
                mg3 = self.get_main_group(group_no_3) #n

                if mg2==mg3:
                    a_nm=0
                else:
                    try:
                # have to do try except as ij is in same row as ji and does not have thei own row
                        try:
                            index=df_ip.index[(df_ip["i"] == mg3) & (df_ip["j"] == mg2)][0]
                            Anm,Bnm,Cnm=df_ip.loc[index,"Aij"],df_ip.loc[index,"Bij"],df_ip.loc[index,"Cij"]

                        except:
                            index=df_ip.index[(df_ip["i"] == mg2) & (df_ip["j"] == mg3)][0]
                            Anm,Bnm,Cnm=df_ip.loc[index,"Aji"],df_ip.loc[index,"Bji"],df_ip.loc[index,"Cji"]
                            
                    except:
                        print(f"Interaction parameters not found for groups {mg2} and {mg3}")
                        return None
                    a_nm=Anm+Bnm*T+Cnm*T**2
                cap_psi_nm=np.exp(-a_nm/(T))
                
                sum_n_1+=cap_theta_n*cap_psi_nm
            
            sum_m_2+=cap_theta_m*cap_psi_km/sum_n_1
        
        ln_kap_gamma_k=cap_q_k*(1-np.log(sum_m_1)-sum_m_2)
        return ln_kap_gamma_k

    #=======================================ln gamma_k^i=========================================
    def ln_cap_gamma_i(self,group_no_1,species_no):
        T=self.T
        list_groups=self.list_groups

        index = df_sv.index[df_sv["No."] == group_no_1][0]
        cap_q_k=df_sv.loc[index,"Q"] #get from group contribution table

        df=list_groups[species_no]

        sum_m_1=0
        sum_m_2=0
        #==========================sum m 1:===========================
        for m in range(len(df)): # sum over all gropues in molecule
            
            group_no_2=df.loc[m,"Group_No"]
            
            cap_theta_m_i=self.cap_theta_i(group_no_2,species_no)

            # convert subgroup number to main group number for interaction parameters
            mg1 = self.get_main_group(group_no_1)
            mg2 = self.get_main_group(group_no_2)

            if mg1==mg2:
                a_mk,a_km=0,0 # self interaction
            else:
                # have to do try except as ij is in same row as ji and does not have thei own row
                try:
                    try:
                        index=df_ip.index[(df_ip["i"] == mg2) & (df_ip["j"] == mg1)][0]
                        Amk,Bmk,Cmk=df_ip.loc[index,"Aij"],df_ip.loc[index,"Bij"],df_ip.loc[index,"Cij"]
                        Akm,Bkm,Ckm=df_ip.loc[index,"Aji"],df_ip.loc[index,"Bji"],df_ip.loc[index,"Cji"]
                    except:
                        index=df_ip.index[(df_ip["i"] == mg1) & (df_ip["j"] == mg2)][0]
                        Amk,Bmk,Cmk=df_ip.loc[index,"Aij"],df_ip.loc[index,"Bij"],df_ip.loc[index,"Cij"]
                        Akm,Bkm,Ckm=df_ip.loc[index,"Aji"],df_ip.loc[index,"Bji"],df_ip.loc[index,"Cji"]
                except:
                    print(f"Interaction parameters not found for groups {mg1} and {mg2}")
                    return None
                a_mk=Amk+Bmk*T+Cmk*T**2
                a_km=Akm+Bkm*T+Ckm*T**2
            cap_psi_mk=np.exp(-a_mk/T)
            cap_psi_km=np.exp(-a_km/T)

            sum_m_1+=cap_theta_m_i*cap_psi_mk

            
            sum_n_1=0
            #==============================sum m 2:=========================
            #sum n:
            for n in range(len(df)): # sum over all grupes in molecule
                group_no_3=df.loc[n,"Group_No"]
                cap_theta_n_i=self.cap_theta_i(group_no_3,species_no)

                # convert subgroup number to main group number for interaction parameters
                mg3 = self.get_main_group(group_no_3)
                if mg2==mg3:
                    a_nm=0
                else:
                    try:
                # have to do try except as ij is in same row as ji and does not have thei own row
                        try:
                            index=df_ip.index[(df_ip["i"] == mg3) & (df_ip["j"] == mg2)][0]
                            Anm,Bnm,Cnm=df_ip.loc[index,"Aij"],df_ip.loc[index,"Bij"],df_ip.loc[index,"Cij"]
                        except:
                            index=df_ip.index[(df_ip["i"] == mg2) & (df_ip["j"] == mg3)][0]
                            Anm,Bnm,Cnm=df_ip.loc[index,"Aji"],df_ip.loc[index,"Bji"],df_ip.loc[index,"Cji"]
                    except:
                        print(f"Interaction parameters not found for groups {mg2} and {mg3}")
                        return None
                    a_nm=Anm+Bnm*T+Cnm*T**2
                cap_psi_nm=np.exp(-a_nm/(T))
                
                sum_n_1+=cap_theta_n_i*cap_psi_nm
            
            sum_m_2+=cap_theta_m_i*cap_psi_km/sum_n_1
        
        ln_kap_gamma_k_i=cap_q_k*(1-np.log(sum_m_1)-sum_m_2)
        return ln_kap_gamma_k_i

    #=======================================residual=============================================
    def ln_gamma_res_i(self,species_no):
        list_groups=self.list_groups
        df=list_groups[species_no]
        
        sum_k=0

        for k in range(len(df)): # sum over all gropupes present in the molecule i
            group_no=df.loc[k,"Group_No"]
            
            index2 = df.index[df["Group_No"] == group_no][0]
            nu_m_i=df.loc[index2, "number"]


            ln_cap_gamma_k_i=self.ln_cap_gamma_i(group_no,species_no)
            ln_cap_gamma_k=self.ln_cap_gamma(group_no)
            sum_k+=nu_m_i*(ln_cap_gamma_k-ln_cap_gamma_k_i)
        return sum_k

    #========================================combinatorial====================================0
    def cap_phi(self,species_no):
        mol_lst=self.mol_lst
        list_groups=self.list_groups
        
        # r_i
        sum_j=0
        def r(species_no):
            sum_k=0
            df2=list_groups[species_no]
            for k in range(len(df2)): #sum over number of groupes in molecule
                #print(f"her1_{j}")
                group_no=df2.loc[k,"Group_No"]
                try:
                    index = df2.index[df2["Group_No"] == group_no][0]
                except:
                    continue
                nu_ik=df2.loc[index, "number"]
                
                index2 = df_sv.index[df_sv["No."] == group_no][0]
                
                cap_r_k=df_sv.loc[index2,"R"]
                
                sum_k+=cap_r_k*nu_ik

            r_i=sum_k
            return r_i
        
        for j in range(len(list_groups)): # sum over number of species in mixture
            
            x_j=mol_lst[j]
            #print(f"her1_{j}")
            r_j=r(j)
            #print(f"her2_{j}")
            sum_j+=r_j*x_j
        
        x_i=mol_lst[species_no]
        r_i=r(species_no)    
        
        cap_phi_i=r_i*x_i/sum_j
        
        return cap_phi_i
    
    def phi_mark(self,species_no):
        mol_lst=self.mol_lst
        list_groups=self.list_groups
        # r_i
        sum_j=0
        def r(species_no):
            sum_k=0
            df2=list_groups[species_no]
            for k in range(len(df2)): #sum over number of groupes in molecule
                
                group_no=df2.loc[k,"Group_No"]
                index = df2.index[df2["Group_No"] == group_no][0]
                nu_ik=df2.loc[index, "number"]

                index2 = df_sv.index[df_sv["No."] == group_no][0]
                cap_r_k=df_sv.loc[index2,"R"]
                sum_k+=cap_r_k*nu_ik

            r_i=sum_k
            return r_i
        
        for j in range(len(list_groups)): # sum over number of species in mixture
            x_j=mol_lst[j]
            r_j=r(j)
            
            sum_j+=r_j**(3/4)*x_j
        
        x_i=mol_lst[species_no]
        r_i=r(species_no)    

        phi_i=r_i**(3/4)*x_i/sum_j
        return phi_i

    def q(self,species_no):
        list_groups=self.list_groups
        df2=list_groups[species_no]
        sum_k=0
        for k in range(len(df2)): #sum over number of groupes in molecule
            
            group_no=df2.loc[k,"Group_No"]
            index = df2.index[df2["Group_No"] == group_no][0]
            nu_ik=df2.loc[index, "number"]

            index2 = df_sv.index[df_sv["No."] == group_no][0]
            cap_q_k=df_sv.loc[index2,"Q"]
            sum_k+=cap_q_k*nu_ik

        q_i=sum_k
        
        return q_i

    def theta(self,species_no):
        list_groups=self.list_groups
        mol_lst=self.mol_lst
        # r_i
        sum_j=0
        
        for j in range(len(list_groups)): # sum over number of species in mixture
            x_j=mol_lst[j]
            q_j=self.q(j)
            
            sum_j+=q_j*x_j
        
        x_i=mol_lst[species_no]
        q_i=self.q(species_no)    

        theta_i=q_i*x_i/sum_j
        return theta_i
    #gamma
    def ln_gamma_comb_i(self,species_no):
        mol_lst=self.mol_lst
        z=10
        #df=list_groups[species_no]
        cap_phi_i = self.cap_phi(species_no)     # volume fraction 
        phi_mark_i= self.phi_mark(species_no)    # modified volume fraction 
        x_i= mol_lst[species_no] # molefraction 
        theta_i= self.theta(species_no)          # surface area fraction 
        q_i=self.q(species_no)

        ln_cap_gamma=np.log(phi_mark_i/x_i)+1-phi_mark_i/x_i-z/2*q_i*(np.log(cap_phi_i/theta_i)+1-cap_phi_i/theta_i)
        return ln_cap_gamma

    #=================================main=================================
    def gamma_singular(self,species_no):
        gamma_i=np.exp(self.ln_gamma_comb_i(species_no)+self.ln_gamma_res_i(species_no))
        return gamma_i
    
    def gamma_total(self):
        mol_lst=self.mol_lst
        lst=[]
        for i in range(len(mol_lst)): 
            gamma_main_i=self.gamma_singular(i)
            lst.append(gamma_main_i)
        return lst





















