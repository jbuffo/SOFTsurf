function SOFTsurf = makeSurface(salinity,slope,C_e,beta,g,k_s,saltname,path)

    % DATA TYPES: salinity - array; slope - double; C_e - double; beta -
    % double; g - double; k_s - double; saltname - char
    
    %% Switching from dimensional to dimensionless parameters for Jamie's code & building input files accordingly
    %% for a range of salinities and domain sizes to produce S_ice vs. dT/dz vs. S_oc surfaces
    
    num_salinities=length(salinity);         % Number of salinities being swept through
    num_sim_size=4;           % Number of domain sizes simulated
    cell_array_hold=cell(num_salinities,num_sim_size);     % preallocate cell array to add results into
    
    %salinity=[130,170];                         % Salinity values to cover
    simulation_size=[10,80,1000,2000];          % y domain size in meters for each domain size run
    
    parfor i=1:num_salinities
    
        system(['mkdir ',num2str(salinity(i)),'ppt']);
    
        %% DOMAIN SIZE & OUTPUT SPECS (Dimensions based on scale height, desired resolution, chk/plt output intervals)
        
        % Some preset values for domain size/resolution sweeps at 10m, 80m, 1km, &
        % 2km domain heights to generate S_ice vs. dT/dz curves
        
        for j=1:num_sim_size
    
            system(['mkdir ',num2str(salinity(i)),'ppt/',num2str(simulation_size(j)),'m']);
            
            if simulation_size(j)==10
                domain_height=1.0;          % In scale height units
                y_grid=256;                 % Number of grid elements in y-direction (resolution=domain_height*h/y_grid)
                x_grid=256;                 % Number of grid elements in x-direction
                chk_interval=10000;          % Checkpoint file interval (number of time steps)
                plot_period=0.05;          % Plot file period (in dimensionless time)
                max_step=25000;             % Max number of time step iterations
                max_time=50.0;              % Max dimensionless time
                min_time=50.0;              % Minimum dimensionless time
                h=10;                       % Scale height
            elseif simulation_size(j)==80
                domain_height=8.0;          % In scale height units
                y_grid=256;                 % Number of grid elements in y-direction (resolution=domain_height*h/y_grid)
                x_grid=32;                  % Number of grid elements in x-direction
                chk_interval=10000;         % Checkpoint file interval (number of time steps)
                plot_period=0.5;           % Plot file period (in dimensionless time)
                max_step=120000;            % Max number of time step iterations
                max_time=50.0;              % Max dimensionless time
                min_time=50.0;              % Minimum dimensionless time
                h=10;                       % Scale height
            elseif simulation_size(j)==1000
                domain_height=1.0;          % In scale height units
                y_grid=128;                 % Number of grid elements in y-direction (resolution=domain_height*h/y_grid)
                x_grid=128;                 % Number of grid elements in x-direction
                chk_interval=100000;         % Checkpoint file interval (number of time steps)
                plot_period=0.005;          % Plot file period (in dimensionless time)
                max_step=500000;            % Max number of time step iterations
                max_time=50.0;              % Max dimensionless time
                min_time=50.0;              % Minimum dimensionless time
                h=1000;                     % Scale height
            elseif simulation_size(j)==2000
                domain_height=2.0;          % In scale height units
                y_grid=128;                 % Number of grid elements in y-direction (resolution=domain_height*h/y_grid)
                x_grid=16;                  % Number of grid elements in x-direction
                chk_interval=100000;         % Checkpoint file interval (number of time steps)
                plot_period=0.005;          % Plot file period (in dimensionless time)
                max_step=850000;            % Max number of time step iterations
                max_time=50.0;              % Max dimensionless time
                min_time=50.0;              % Minimum dimensionless time
                h=1000;                     % Scale height
            else
            end
            
            %% INPUTS (Dimensional parameters for the problem)
            
            T_s=-173.15;                    % Surface temperature (C)
            C_i=salinity(i);                   % Initial ocean composition (ppt)
            %slope=-0.0913;                  % Linear liquidus slope (-0.0913 for NaCl, -0.022857 for MgSO4)
            T_oc=0+slope*C_i+0.01;          % Ocean temperature (C)
            %C_e=230;                        % Eutectic composition (ppt) (230 for NaCl, 175 for MgSO4)
            T_e=0+slope*C_e;                % Eutectic temperature (C)
            phi_s=0;                        % Top BC porosity
            phi_oc=1;                       % Bottom BC porosity
            c_i=2000;                       % Specific heat of ice (J/kg*K)
            c_br=3985;                      % Specific heat of ocean (J/kg*K)
            rho_i=917;                      % Density of ice (kg/m^3)
            %beta=5.836*10^-4;%7.7*10^-4;    % Density coefficient for salt (kg/g) (7.7*10^-4 NaCl with range 6.4 to 7.8)
                                            % (8.3*10^-4 MgSO4 with range 6.6-10) from Zhu et al 2017
            rho_br=1000+1000*beta*C_i;      % Density of ocean (kg/m^3)
            k_i=2.0;                        % Thermal conductivity of ice (W/m*K)
            k_br=0.6;                       % Thermal conductivity of ocean (W/m*K)
            %h=10;                          % Scale height (m)
            L=334774;                       % Latent heat of fusion (J/kg) (334774 for ice-water transition)
            eta=1.88*10^-3;                 % Dynamic viscosity (Pa*s) (1.88*10^-3 for water near 0C)
            %g=1.32;                         % Gravity (m/s^2) (1.32 for Europa) (0.113 for Enceladus)
            %k_s=2*10^-9;                    % Salt diffusivity in water (m^2/s) (2*10^-9 for typical runs)
            p_s=0.001;                      % Partition coefficient
            K_0=20*10^-10;                  % Permeability factor (m^2)
            d=5*10^-5;                      % Hele-Shaw cell spacing (m)
            alpha=2.10*10^-4;               % Thermal expansion coefficient for water
            
            %% CALCULATIONS/OUTPUTS
            
            del_T=(0+slope*C_i)-T_e;                            % Characteristic temperature scale
            del_C=C_e-C_i;                                      % Charactersitic composition scale
            thetaT_top=(T_s-T_e)/del_T;                         % Dimensionless temperature @ the surface
            thetaT_oc=(T_oc-T_e)/del_T;                         % Dimensionless temperature in the ocean
            thetaC_top=(0-C_e)/del_C;                           % Dimensionless composition @ the surface
            thetaC_oc=(C_i-C_e)/del_C;                          % Dimensionless composition in the ocean
            St=L/(c_br*del_T);                                  % Stefan number
            c_p=c_i/c_br;                                       % Specific heat ratio
            H_top=St*phi_s+(phi_s+(1-phi_s)*c_p)*thetaT_top;    % Enthalpy @ surface
            H_oc=St*phi_oc+(phi_oc+(1-phi_oc)*c_p)*thetaT_oc;   % Enthalpy in ocean
            kappa_i=k_i/(rho_i*c_i);                            % Thermal diffusivity of ice
            kappa_br=k_br/(rho_br*c_br);                        % Thermal diffusivity of ocean
            Pr=eta/(rho_br*kappa_br);                           % Prandtl number
            Ra_C=(rho_br*g*beta*del_C*h^3)/(kappa_br*eta);      % Rayleigh number for composition
            Ra_T=(rho_br*g*alpha*del_T*h^3)/(kappa_br*eta);     % Rayleigh number for temperature
            Tau=(h^2)/kappa_br;                                 % Dimensionless timescale (diffusive)
            K_ratio=kappa_i/kappa_br;                           % Ratio of heat diffusivities
            C_ratio=((1-p_s)*C_e)/del_C;                        % Composition ratio
            k_ratio=k_i/k_br;                                   % Conductivity ratio
            Le=kappa_br/k_s;                                    % Lewis number
            R=(12*K_0)/d^2;                                     % Reluctance
            Da=K_0/h^2;                                         % Darcy number
            Rm_C=Da*Ra_C;                                       % Mushy Layer number (composition)
            Rm_T=Da*Ra_T;                                       % Mushy Layer number (temperature)
            
            %% DISPLAY VALUES
            
            % Vars={'thetaT_top', 'thetaT_oc', 'thetaC_oc', 'St', 'c_p', 'H_top', 'H_oc', 'Tau', 'K_ratio', 'C_ratio', 'k_ratio',...
            %      'Le', 'R', 'Da', 'Rm_C', 'Rm_T', 'T_e'}';
            % Vals=[thetaT_top, thetaT_oc, thetaC_oc, St, c_p, H_top, H_oc, Tau, K_ratio, C_ratio, k_ratio,...
            %      Le, R, Da, Rm_C, Rm_T, T_e]';
            % Out=table(Vars,Vals);
            
            %% MODIFY INPUT FILE
            
            func_replace_string('inputs_0',['inputs',num2str(i)],'bc.bulkConcentrationLoVal=',['bc.bulkConcentrationLoVal=0.0 ',num2str(thetaC_oc)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'bc.enthalpyHiVal=',['bc.enthalpyHiVal=0.0 ',num2str(H_top)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'bc.enthalpyLoVal=',['bc.enthalpyLoVal=0.0 ',num2str(H_oc)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'bc.liquidConcentrationLoVal=',['bc.liquidConcentrationLoVal=0.0 ',num2str(thetaC_oc)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'bc.temperatureHiVal=',['bc.temperatureHiVal=0.0 ',num2str(thetaT_top)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'bc.temperatureLoVal=',['bc.temperatureLoVal=0.0 ',num2str(thetaT_oc)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'main.checkpoint_interval=',['main.checkpoint_interval=',num2str(chk_interval)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'main.domain_height=',['main.domain_height=',num2str(domain_height)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'main.max_step=',['main.max_step=',num2str(max_step)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'main.max_time=',['main.max_time=',num2str(max_time)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'main.min_time=',['main.min_time=',num2str(min_time)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'main.num_cells=',['main.num_cells=',num2str(x_grid),' ',num2str(y_grid),' 32'])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'main.plot_period=',['main.plot_period=',num2str(plot_period)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.K=',['parameters.K=',num2str(K_ratio)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.compositionRatio=',['parameters.compositionRatio=',num2str(C_ratio)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.eutecticComposition=',['parameters.eutecticComposition=',num2str(C_e)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.eutecticTemp=',['parameters.eutecticTemp=',num2str(T_e)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.heatConductivityRatio=',['parameters.heatConductivityRatio=',num2str(k_ratio)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.initialComposition=',['parameters.initialComposition=',num2str(C_i)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.lewis=',['parameters.lewis=',num2str(Le)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.liquidusSlope=',['parameters.liquidusSlope=',num2str(slope)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.nonDimReluctance=',['parameters.nonDimReluctance=',num2str(R)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.rayleighComp=',['parameters.rayleighComp=',num2str(Rm_C)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.rayleighTemp=',['parameters.rayleighTemp=',num2str(Rm_T)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.referenceSalinity=',['parameters.referenceSalinity=',num2str(C_e)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.referenceTemperature=',['parameters.referenceTemperature=',num2str(T_e)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.specificHeatRatio=',['parameters.specificHeatRatio=',num2str(c_p)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'parameters.stefan=',['parameters.stefan=',num2str(St)])
            func_replace_string(['inputs',num2str(i)],['inputs',num2str(i)],'main.output_folder=',['main.output_folder=./',num2str(salinity(i)),'ppt/',num2str(simulation_size(j)),'m'])
            
            %% MAKE INPUT FILE AN EXECUTABLE AND MOVE TO ASSOCIATED SALINITY AND DOMAIN SIZE REPO.
            
            system(['chmod +x inputs',num2str(i)]);
            system(['mv inputs',num2str(i),' ',num2str(salinity(i)),'ppt/',num2str(simulation_size(j)),'m']);
    
            %% RUN SOFTBALL SIMULATION
            system(['LD_LIBRARY_PATH=/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH ',path,'/mushy-layer/execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex ./',num2str(salinity(i)),'ppt/',num2str(simulation_size(j)),'m/inputs',num2str(i)]);
    
            %% ANALYZE LAST HDF5 FILE TO EXTRACT S_ICE VS. S_OC VS. DT/DZ
            % grab stuff using Jamie's code
            addpath(genpath([path,'/mushy-layer/matlab']));
            output_dir = ['./',num2str(salinity(i)),'ppt/',num2str(simulation_size(j)),'m'];
            actual_plot_prefix = 'plt';
            thisFrame = max_step;
            ml = MushyLayerOutput(2, thisFrame, output_dir, actual_plot_prefix, true);
    
            % extract values and dimensionalize
            S  = ml.dataForComp(ml.components.Bulkconcentration).';
            porosity = ml.dataForComp(ml.components.Porosity).';
            dim_S = flipud(S*(C_e-C_i)+C_e);
            Phi = flipud(porosity);
            avg_dim_S = mean(dim_S,2);
            avg_Phi = mean(Phi,2);
    
            % find cells above mushy layer
                % preset depth_interface in case run failed and outputs are
                % garbage to avoid errors
            depth_interface=1;
            for k=1:length(avg_dim_S)
                if avg_Phi(k)<10^-10 && k<0.85*length(avg_dim_S)
                    depth_interface=k;
                else
                    break
                end
            end
    
            % calculate depths and thermal gradients
            resolution=simulation_size(j)/y_grid;
            depths=[resolution/2:resolution:depth_interface*resolution-resolution/2]';
            dT_dz=(T_oc-T_s)./depths;
    
            % save values to cell array {S_ocean, dT/dz, S_ice}
            S_ocean=C_i*ones(depth_interface,1);
            dT_dz_save=dT_dz(1:depth_interface);
            S_ice_save=avg_dim_S(1:depth_interface);
    
            cell_array_hold{i,j}=[S_ocean,dT_dz_save,S_ice_save];
    
            % delete folder for that salinity-domain size pair as we're done with it and can save space
            system(['rm -rf ',num2str(salinity(i)),'ppt/',num2str(simulation_size(j)),'m']);
        
        end
    
    end
    
    %% OUTSIDE OF PARFOR LOOP - SAVE ALL DATA INTO CELL ARRAY
    save('All_values_array.mat','cell_array_hold');
    
    %% ORGANIZE DATA FOR CURVE FITTING AND INTERPOLATION
    % matrices for gathering fitlines and scatter points for plotting and
    % interpolation
    fitlines=cell(1,num_salinities);
    fitpoints=cell(1,num_salinities);
    actualpoints=cell(1,num_salinities);
    
    for i=1:num_salinities
        hold_mat=[];
        for j=1:num_sim_size
            hold_mat=[hold_mat;cell_array_hold{i,j}];
        end
    
        % hold_mat == [S_oc, dT/dz, S_ice]
        sort_hold_mat=sortrows(hold_mat,2);
    
        % curve fit to quasi-analytical solution from Buffo et al., 2021 paper
        % fit from 2021 paper
        ft=fittype('a+(b*(dT+c)/(d+f*dT))*(1-h*exp(-k/dT))','dependent',{'y'},'independent',...
            {'dT'},'coefficients',{'a','b','c','d','f','h','k'});

        % bounds for coefficients to prevent fit errors
        opts = fitoptions('Method', 'NonlinearLeastSquares');
        opts.Lower = [-Inf, -Inf, -Inf, 1e-5, -Inf, -Inf, 0];  % avoid d near 0, k ≥ 0
        opts.Upper = [Inf, Inf, Inf, Inf, Inf, Inf, Inf];

        % turn fit warnings off
        warning('off', 'all');

        % fit curve (250x) and create values along the best fit line (max r^2)
        yfit_hold=[];
        rsq_hold=[];
        f_hold=cell(1,250);
        for k=1:250
            f=fit(sort_hold_mat(:,2),sort_hold_mat(:,3),ft,opts);
            xfit = sort([linspace(0,1,10000), linspace(min(sort_hold_mat(:,2)),max(sort_hold_mat(:,2)),100000)]);        % add 0-1 values
            yfit = feval(f, xfit);                                                              % prevent negative salt(?)
            rsq = 1 - sum((sort_hold_mat(:,3) - f(sort_hold_mat(:,2))).^2) / sum((sort_hold_mat(:,3) - mean(sort_hold_mat(:,3))).^2);
            yfit_hold=[yfit_hold yfit];
            rsq_hold=[rsq_hold rsq];
            f_hold{k}=f;
        end
        [maxval, loc]=max(rsq_hold);
        yfit=yfit_hold(:,loc);
        f=f_hold{loc};
        % force monotonicity with increasing ocean salinity outside of real data
        % (remedy decreasing/negative fit values for high S_oc & low dT/dz beyond lowest simulated values)
        for l=1:length(yfit)
            if i==1
                continue
            elseif xfit(l)>=min(sort_hold_mat(:,2))
                continue
            else
                f_test=[];
                for xx=1:i-1
                    f_test_fct=fitlines{xx};
                    f_test=[f_test feval(f_test_fct,xfit(l))];
                end
                if yfit(l)<max(f_test)
                    yfit(l)=max(f_test);
                else
                    continue
                end
            end
        end
        % append values and fit lines to cell arrays for plotting and
        % interpolation
        fitlines{i}=f;
        fitpoints{i}=[salinity(i)*ones(length(xfit),1),xfit',yfit];
        actualpoints{i}=sort_hold_mat;
    end
    
    %% CREATE INTERPOLATED SURFACE OF THE RESULTANT FIT LINES (WITH MAX S_ICE=S_OC) AND PLOT TO VALIDATE
    % plot results to make sure the fits/interpolation looks good
    h = figure;
    hold on
    all_vals_fit=[];
    all_vals_real=[];
    % plot all the fit lines for each S_oc value simulated and accumulate all
    % scatter values for the interpolation
    for i=1:num_salinities
        hold_mat_fit=fitpoints{i};
        hold_mat_real=actualpoints{i};
        plot3(hold_mat_fit(:,1),hold_mat_fit(:,2),hold_mat_fit(:,3),'r','LineWidth',2)
        all_vals_fit=[all_vals_fit;hold_mat_fit];
        all_vals_real=[all_vals_real;hold_mat_real];
    end
    % plot all simulated datapoints
    scatter3(all_vals_real(:,1),all_vals_real(:,2),all_vals_real(:,3),'k')
    % create mesh grid to sample interpolation surface on

    %[xq,yq]=meshgrid([0:5:max(salinity)],[min(all_vals_fit(:,2)):5:max(all_vals_fit(:,2))]);
    [xq,yq]=meshgrid([0:2.5:C_e],[0:2.5:max(all_vals_fit(:,2))]);

    dT_vals=[0:2.5:max(all_vals_real(:,2))]';
    % add S_oc=0 => S_ice=0 values
    all_vals_fit=[all_vals_fit;zeros(length(dT_vals),1),dT_vals,zeros(length(dT_vals),1)];
    % add S_oc=C_e => S_ice=C_e values
    all_vals_fit=[all_vals_fit;C_e*ones(length(dT_vals),1),dT_vals,C_e*ones(length(dT_vals),1)];
    % interpolate gridded data points and restrict S_ice<S_oc
    zq=griddata(all_vals_fit(:,1),all_vals_fit(:,2),min(all_vals_fit(:,3),all_vals_fit(:,1)),xq,yq,'linear');
    surf(xq,yq,zq);
    xlim([0 C_e])
    ylim([0 200])
    zlim([0 C_e+10])
    xlabel('Ocean Salinity (ppt)')
    ylabel('Thermal Gradient (K/m)')
    zlabel('Ice Salinity (ppt)')
    set(gca,'FontSize',20)
    
    % save the figure
    savefig(h,'Surface.fig')
    
    %% NOW SAVE THIS INTERPOLATED SURFACE AS A FUNCTION HANDLE THAT CAN BE CALLED AS S_ICE=GETVAL(S_OC,DT/TZ)
    SOFTsurf.S_ice = @(S_oc,dTdz) griddata(all_vals_fit(:,1),all_vals_fit(:,2),min(all_vals_fit(:,3),all_vals_fit(:,1)),S_oc,dTdz);

    % saving a gridded interpolant as an option to speed up search (test % errors more thoroughly - so far < 1.5% for all tested values)
    [xq2,yq2]=meshgrid([0:2.5:C_e],[0:0.1:max(all_vals_fit(:,2))]);
    zq2=griddata(all_vals_fit(:,1),all_vals_fit(:,2),min(all_vals_fit(:,3),all_vals_fit(:,1)),xq2,yq2,'linear');
    interp_grid = griddedInterpolant(xq2',yq2',zq2','linear');
    SOFTsurf.S_ice_fast = @(S_oc,dTdz) interp_grid(S_oc,dTdz);
    
    %% SAVE THE SURFACE AS A .MAT SO THAT IT CAN BE RELOADED AND UTILIZED
    save(['Surf',saltname,'.mat'],'SOFTsurf','salinity','slope','C_e','beta','g','k_s','saltname');

    %% ADD COPIES OF ALL SAVED FILES TO A NEW DIRECTORY SO IT CAN BE UPLOADED TO GITHUB, BUILDING A LIBRARY OF SALINE ICE CHEMISTRIES FOR DIFFERENT BINARY SALTS
    system(['mkdir ',saltname]);
    system(['mv All_values_array.mat ',saltname,'/']);
    system(['mv Surface.fig ',saltname,'/']);
    system(['mv Surf',saltname,'.mat ',saltname,'/']);
    system(['mv ',saltname,' Salts/']);

    %% CLEAN UP THE PARENT REPO - FILES WILL BE IN 'saltname' FOLDER
    for i=1:num_salinities
        system(['rm -rf ',num2str(salinity(i)),'ppt/'])
    end
    system('rm -f diagnostics.csv')
    system('rm -f diagnosticsLatest.csv')
    system('rm -f pout.0')
    
end
