#include "params.cpp"
#include "discequations.cpp"
#include "postproc.cpp"
#include "allocate_arrays.cpp"





int main()
{

double residual=1;

int count=0;
bool loop_switch=true;

init_variables(u,v,p);

while(loop_switch)
{

 	solve_x_mom(u,unp1,v,p);
 	solve_y_mom(v,vnp1,u,p);
 	solve_pressure(u,v,p,pnp1);
 	apply_bc(u,v,p);
	apply_bc(unp1,vnp1,pnp1);
	residual=compute_residual(unp1,vnp1);
	swap_variables(u,v,p,unp1,vnp1,pnp1);
	
	
	if(residual<tol&&(count!=0))	
	{
	compute_collocated_values(unp1,vnp1,pnp1,uc,vc,pc);
	write_file(uc,vc,pc,nx,ny);
	std::cout<<"converged to a tolerance of: "<<tol<<" in "	<<count<<" Iterations"<<std::endl;
	loop_switch=false;
	
	}
	
	if(count%print_interval==0)
	{
	compute_collocated_values(unp1,vnp1,pnp1,uc,vc,pc);
	write_file(uc,vc,pc,nx,ny);
	std::cout<<"|| iteration is: "<<count<<", continuity residual : " <<residual<<" ||"<<std::endl;
	
	}
	
	
	count+=1;	

}

return 0;

}
