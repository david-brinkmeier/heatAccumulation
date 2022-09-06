classdef materialprops
    % Generate structure to store relevant material constants for a material
    % types of properties: https://www.mathworks.com/help/matlab/matlab_oop/mutable-and-immutable-properties.html
    properties
        % thermophysical quantities
        name            string % Name as string
        lambda          double % Heat conductivity in W/(m*K) | kappa = lambda_therm/(rho*cp)
        rho             double % Density in kg/m³
        cp              double % Specific heat in J/kgK
        Tcritical       double % Critical temperature in °C
        eta_abs         double % Residual heat (eta_abs*Ep = energy input) ||| | Weber 0.125 / Faas 0.55*0.38 = 0.209
        hth             double % Threshold for the absorbed fluence in J/cm² ||| (definition as in Holder 2021 - 10.1007/s00339-021-04455-3)
        hvL             double % J/mm^3 Volume specific melt energy
        hV              double % J/mm^3 Volume specific evoparization energy: hth = lalpha*10^-6*hV*10^2; % J/cm²
        % optical quantities
        wavelength      double % Associated optical wavelength
        refrindex(1,2)  cell   % [n,k]. If k is not provided, it is derived from linear absorption coefficient and optical wavelength (see docs)
        pdepth          double = [] % optical penetration depth in m, if not explicitly provided calculated as 1/linabscoeff
        linabscoeff     double = [] % linear absorption coefficient; can be derived from extinction coefficient (see docs)
        A0              double = [] % Absorption coefficient at normal incidence. If not provided it is internally calculated from optical properties.
    end
    properties (SetAccess = protected)
        lthermal    double % Thermal diffusion length in m [Graf Laser 2, pp. 158)
    end
    properties (Dependent) %(Dependent, SetAccess = private)
        kappa       double % Heat diffusivity in m²/s
        hth_hV      double % Threshold fluence in J/cm^2 calculated from other material properties
    end
    
    methods
        % Constructor
        function obj = materialprops(varargin)
            if ~isempty(varargin) && numel(varargin) == 1
                if strcmpi('crni_1030nm',varargin{:})
                    % thermophysical / laser processing quantities
                    obj.name = varargin{:};         % set Name
                    obj.lambda = 15;                % Heat conductivity in W/(m*K)
                    obj.rho = 7900;                 % Density in kg/m³
                    obj.cp = 477;                   % Specific heat in J/(kg*K)
                    obj.Tcritical = 1500-20;        % Critical delta_T
                    obj.eta_abs = 0.125;            % Residual heat (eta_abs*Ep = energy input)
                    obj.hth = 0.12;                 % Threshold for the absorbed fluence in J/cm²
                    obj.hvL = 10;                   % J/mm^3 melt energy
                    obj.hV = 61;                    % J/mm^3 evaporization energy
                    % optical quantities, src https://refractiveindex.info/?shelf=3d&book=metals&page=iron
                    obj.wavelength = 1030e-9;       % associated optical radiation wavelength
                    obj.A0 = [];                    % Absorption coefficient at normal indidence. Determined using Fresnel. Set value to override.
                    obj.refrindex = {2.9421,[]};    % n+ik; k may be calculated from other properties
                    obj.pdepth = 20e-9;             % Energy penetration depth in m
													% pdepth src: https://doi.org/10.1007/s00339-021-04455-3
                    obj.linabscoeff = 4.7696e+07;   % optical linear absorption coefficient in SI units - 1/m
                end
                if strcmpi('bk7_1030nm',varargin{:})
                    obj.name = varargin{:};         % set Name
                    obj.lambda = 1.1;               % Heat conductivity in W/(m*K)
                    obj.rho = 2430;                 % Density in kg/m³
                    obj.cp = 840;                   % Specific heat in J/(kg*K)
                    obj.Tcritical = 884-22;         % Critical delta_T
                    obj.eta_abs = 0.5;              % Residual heat (eta_abs*Ep = energy input)
                    obj.hth = NaN;                  % Threshold fluence in J/cm²
                    obj.hvL = NaN;                  % J/mm^3 melt energy
                    obj.hV = NaN;                   % J/mm^3 evaporization energy
                    obj.pdepth = [];                % Energy penetration depth in nm
                    % optical quantities, src https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
                    obj.wavelength = 1030e-9;       % associated optical radiation wavelength
                    obj.A0 = [];                    % Absorption coefficient at normal indidence. Determined using Fresnel. Set value to override.
                    obj.refrindex = {1.5071,[]};    % n+ik; k may be calculated from other properties
                    obj.pdepth = [];                % Energy penetration depth in m
                    obj.linabscoeff = 0.12245;      % optical linear absorption coefficient in SI units - 1/m
                end
            else
                obj.name = 'unknown';
            end
        end
        
        %% Calculate protected properties

        % Calculate thermal diffusion length lthermal which depends on
        % material property kappa and time for diffusion (for example 1/frep)
        function obj = get_lthermal(obj,timestep)
            % usage direct access: get_lthermal(mat,laser.tpp).lthermal
			% distance where local temperature reaches ~9% of initial surface temperature
            obj.lthermal = 2*sqrt(obj.kappa*timestep);
        end
        
        %% Set dependent properties using set/get
        % calculate values of dependent properties
        function value = get.kappa(obj)
            value = obj.lambda/(obj.rho*obj.cp); % unit: m^2/s
        end
        
        function value = get.pdepth(obj)
            if isempty(obj.pdepth) && ~isempty(obj.linabscoeff)
                value = obj.linabscoeff^-1;
            else
                value = obj.pdepth;
            end
        end
        
        function value = get.A0(obj)
            % Fresnel: huegel/graf laser in der fertigung (3.22) (3.23) (3.24)
            if isempty(obj.A0)
                try
                    nik_1 = 1+0i; % refractive index air
                    nik_2 = obj.refrindex{1}+1i*obj.refrindex{2}; % refractive index of material
                    r_0 = (nik_1-nik_2)/(nik_1+nik_2);            % amplitude reflection coefficient for 0°
                    R0 = r_0*conj(r_0);                           % reflectivity
                    value = 1-R0;                                 % Absortion = 1-R
                catch
                    warning('Either provide A0 directly OR provide means to calculate complex refractive index of material (-> wavelength and linabscoeff!).')
                    value = [];
                end
            else
                value = obj.A0;
            end
        end
        
        function value = get.refrindex(obj)
            value = obj.refrindex;
            if isempty(value{2}) && ~isempty(obj.linabscoeff) && ~isempty(obj.wavelength)
                value{2} = obj.linabscoeff*obj.wavelength/(4*pi);
            elseif isempty(value{2})
                warning('Complex refractive index "k": Linear absorption coefficient and wavelength required.')
            end
        end
        
        function obj = set.kappa(obj,~)
            warning('kappa property is automatically calculated from other material properties \n');
            fprintf('%s%d\n','thermal diffusivity in m^2/s is: ',obj.kappa)
        end
        
        function value = get.hth_hV(obj)
            if ~isempty(obj.pdepth)
                % if user supplied penetration depth it takes precedence
                value = obj.pdepth*1e3*obj.hV*1e2;
            elseif ~isempty(obj.linabscoeff)
                % else calculate from linear absorption coefficient
                value = (obj.linabscoeff^-1)*1e3*mat.hV*1e2;
            else
                warning('hth calculated from hV is only defined if penetration depth or linear absorption coefficient is provided.')
            end
        end
        
        function fresnel(obj)
            if ~any(cellfun(@isempty,obj.refrindex))
                % SOURCE: huegel/graf laser in der fertigung (3.22) (3.23) (3.24)
                nik_1 = 1+0i; % refractive index air
                nik_2 = obj.refrindex{1}+1i*obj.refrindex{2}; % refractive index of material
                % normal indidence
                r_0 = (nik_1-nik_2)/(nik_1+nik_2); % amplitude reflection coefficient normal incidence
                R0 = r_0*conj(r_0); % reflectivitynormal incidence
                % / normal indidence
                len = 90; % eval points
                theta = linspace(pi/2,0,len); % incidence angles
                % calc fresnel
                rs = -(cos(theta)-sqrt(nik_2^2-(sin(theta)).^2)).^2/...
                    (nik_2^2-1); % reflection coefficient for s-pol
                rp = (nik_2^2.*cos(theta)-sqrt(nik_2^2-(sin(theta)).^2))./...
                    (nik_2^2.*cos(theta)+sqrt(nik_2^2-(sin(theta)).^2)); % reflection coefficient for p-pol
                Rs = rs.*conj(rs); % reflectivity for s-pol, equal to abs(rs^2)
                Rp = rp.*conj(rp); % reflectivity for p-pol, , equal to abs(rp^2)
                Rc = (Rs+Rp)/2; % reflectivity for c-pol
                max = rad2deg(atan(abs((nik_2)))); % (max absorption / brewster)
                % / calc fresnel
                angleofincidence = rad2deg(linspace(pi/2,0,len)); % vector for plotting; A
                fig = genORselectfigbyname(sprintf('Fresnel: %s',obj.name)); clf(fig); ax = axes(fig); widths = [1.5,1];
                h(1) = plot(ax,angleofincidence,1-Rp,'--k','DisplayName','A_{p}','LineWidth',widths(1)); hold on
                h(2) = plot(ax,angleofincidence,1-Rc,'-k','DisplayName','A_{c}','LineWidth',widths(1));
                h(3) = plot(ax,angleofincidence,1-Rs,'-.k','DisplayName','A_{s}','LineWidth',widths(1));
                h(4) = plot(ax,[max,max],[0,1],'--r','LineWidth',1);
                text(ax,max+1,0.95,[num2str(max,'%2.1f'),'°'],'Color','red','FontSize',12)
                text(1,(1-R0)+4/90,sprintf('A_{0} = %1.2f',1-R0),'FontSize',12)
                legend(h(1:3),'Location','Northwest','LineWidth',widths(2));
                title(ax,sprintf('%s',obj.name),'Interpreter','none')
                xlabel(ax,'angle of incidence in °'), ylabel(ax,'absorption coefficient')
                ax.FontSize = 12; ax.LineWidth = 1;
                xlim(ax,[0 90])
                drawnow
            else
                warning('Fresnel plot requires complex refractive index and / or extinction coefficient.')
            end
        end
        
    end
    
    methods (Static)
    end
    
end