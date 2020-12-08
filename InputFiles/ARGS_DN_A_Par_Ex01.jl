# Kinetic parameters
p = Dict([
    :g  => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
	:xE => 0.0607,    # GEM iSynTF concentration (nM)
	:hE => 1.0,       # E2 concentration (nM)
	:mE => 0.0159,    # Maximum synthesis rate regulated by GEM
	:kXE => 139.0,    # GEM:E2 dissociation constant
	:bE => 1.02e-4,   # Fraction of the inactive GEM in the nucleus
	:aE => 0.0048,    # Basal expression (in the absence of GEM)
	:nE => 1.4,       # Hill coefficient for GEM regulation
	:kE => 0.0073,    # GEM:promoter dissociation constant
	:xP => 0.3042,    # Z3PM iSynTF concentration (nM)
	:hP => 1.0,       # E2 concentration (nM)
	:mP => 3.1445,    # Maximum synthesis rate regulated by Z3PM
	:kXP => 32.5,     # Z3PM:Pg dissociation constant
	:bP => 0.0317,    # Fraction of the inactive Z3PM in the nucleus
	:aP => 0.0153,    # Basal expression (in the absence of Z3PM)
	:nP => 2.64,      # Hill coefficient for Z3PM regulation
	:kP => 0.0683,    # Z3PM:promoter dissociation constant
	:mA => NaN,	      # A synthesis rate (nM/min)
	:mD => NaN,	      # D synthesis rate (nM/min)
	:mN => NaN,	      # N synthesis rate (nM/min)
	:mY => 0.1,       # Y maximum synthesis rate (nM/min)
	:a  => 0.01,      # Basal activity
	:n  => 2.5,       # Hill function
	:k  => 1.0,       # Dissociation constant
	:eP => 0.0375,    # A:D association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
	:eM => 0.0001,    # A:D dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
	:bP => 0.375,     # A:N association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
	:bM => 0.0001,    # A:N dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
]);

# Calculate/update synthesis rate according to gene regulation (e.g. iSynTF_mu(x,h,kX,b,m,a,n,k)):
function pSynth(p,iSynTF_mu)
	p[:mA] = iSynTF_mu(p[:xE],p[:hE],p[:kXE],p[:bE],p[:mE],p[:aE],p[:nE],p[:kE]);
	p[:mD] = 0.1;
	p[:mN] = iSynTF_mu(p[:xP],p[:hP],p[:kXP],p[:bP],p[:mP],p[:aP],p[:nP],p[:kP]);
	return;
end