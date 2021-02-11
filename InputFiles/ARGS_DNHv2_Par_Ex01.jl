# Kinetic parameters
p = Dict([
    :g  => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
	:xE => 0.0607,    # GEM iSynTF concentration (nM)
	:hE => 36.0,      # E2 concentration (nM)
	:mE => 0.0159,    # Maximum synthesis rate regulated by GEM
	:kXE => 139.0,    # GEM:E2 dissociation constant
	:bE => 1.02e-4,   # Fraction of the inactive GEM in the nucleus
	:aE => 0.0048,    # Basal expression (in the absence of GEM)
	:nE => 1.4,       # Hill coefficient for GEM regulation
	:kE => 0.0073,    # GEM:promoter dissociation constant
	:xP => 0.3042,    # Z3PM iSynTF concentration (nM)
	:hP => 1.0,       # E2 concentration (nM)
	:mP => 0.00944,   # Maximum synthesis rate regulated by Z3PM
	:kXP => 32.5,     # Z3PM:Pg dissociation constant
	:bP => 0.0317,    # Fraction of the inactive Z3PM in the nucleus
	:aP => 0.0153,    # Basal expression (in the absence of Z3PM)
	:nP => 2.64,      # Hill coefficient for Z3PM regulation
	:kP => 0.0683,    # Z3PM:promoter dissociation constant
	:mY => 0.039917,  # Y maximum synthesis rate (nM/min)
	:aY => 0.0082908, # Basal activity
	:nY => 1.0,       # Hill function
	:kY => 0.018864,  # Dissociation constant
    :bM => 0.01,      # A:N or D:N dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
	:eM => 0.01,      # A:D dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
	:bC1 => 0.05,    # A:N association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
	:bC2 => 0.05,    # D:N association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
	:eC1 => 0.005,    # A:D association (binding) rate - 154B_37A
	:eC2 => 0.005,    # A:D association (binding) rate - 154B_37B
	:eC3 => 0.005,    # A:D association (binding) rate - 154B_155A
	:eC4 => 0.005,    # A:D association (binding) rate - 37A:154B
	:eC5 => 0.005,    # A:D association (binding) rate - 13B:154B
	:eC6 => 0.005,    # A:D association (binding) rate - 37B:154B
	:mA => NaN,	      # BY FUNCTION: A synthesis rate (nM/min)
	:mD => NaN,	      # BY FUNCTION: D synthesis rate (nM/min)
	:mN => NaN,	      # BY FUNCTION: N synthesis rate (nM/min)
	:bP => NaN,       # BY RULE: A:N or D:N association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
	:eP => NaN,       # BY RULE: A:D association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
]);

# Calculate/update synthesis rate according to gene regulation (e.g. iSynTF_mu(x,h,kX,b,m,a,n,k)):
function pSynth(p,iSynTF_mu)
	p[:mA] = iSynTF_mu(p[:xE],p[:hE],p[:kXE],p[:bE],p[:mE],p[:aE],p[:nE],p[:kE]);
	p[:mD] = iSynTF_mu(p[:xE],p[:hE],p[:kXE],p[:bE],p[:mE],p[:aE],p[:nE],p[:kE]);
	p[:mN] = iSynTF_mu(p[:xP],p[:hP],p[:kXP],p[:bP],p[:mP],p[:aP],p[:nP],p[:kP]);
	return;
end
