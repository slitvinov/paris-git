#! /bin/awk -f

BEGIN {abs=1e10; paris_vol = 0e0; PI=3.14159265358979; tolerance = 1e-3}
	{($1 - 1000*1.0e-4) < O ? abs = (1000*1.0e-4 - $1) : abs = ($1 - 1000*1.0e-4); 
	if (abs<1e-5) { paris_vol = $11;
		vol_theory = 1000*1.0e-4*0.05*4+PI*0.2*0.2;
		printf("Vol Paris: %10.6f \n",paris_vol);
		printf("Vol theory: %10.6f \n",vol_theory);
		if (paris_vol > vol_theory) {if ((paris_vol-vol_theory)/vol_theory < tolerance) print "\033[1;32mPASS\033[0m";
							else print "\033[1;31mFAIL\033[0m"}
			else {if ((vol_theory-paris_vol)/vol_theory < tolerance) print "\033[1;32mPASS\033[0m";
							else print "\033[1;31mFAIL\033[0m"}
		}					
	}
END {}
