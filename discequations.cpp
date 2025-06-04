void init_variables(dmatrix &u,dmatrix&v,dmatrix &p);
void solve_x_mom(dmatrix&u,dmatrix&unp1,dmatrix	&v,dmatrix &p);
void solve_y_mom(dmatrix&v,dmatrix &vnp1,dmatrix &u,dmatrix &p);
void solve_pressure(dmatrix&u,dmatrix &v,dmatrix &p,dmatrix&pnp1);
void apply_bc(dmatrix&u,dmatrix &v,dmatrix &p);
double compute_residual(dmatrix&u,dmatrix &v);
void compute_collocated_values(dmatrix &u,dmatrix &v ,dmatrix &p,dmatrix &uc,dmatrix &vc ,dmatrix &pc);
void swap_variables(dmatrix &u,dmatrix &v ,dmatrix &p,dmatrix &unp1,dmatrix &vnp1,dmatrix &pnp1);






void init_variables(dmatrix &u,dmatrix&v,dmatrix &p)
{


for(int i=0;i<=(nx-1);i++)
{
	for(int j=0;j<=ny;j++)
	{
		u(i,j)=0.0;
		u(i,ny)=1.0;	
		u(i,ny-1)=1.0;
		
	}
}
	
for(int i=0;i<=nx;i++)
{
	for(int j=0;j<=(ny-1);j++)
	{	
		v(i,j)=0.0;
		
	}
}


for(int i=0;i<=nx;i++)
{
	for(int j=0;j<=ny;j++)
	{	
		p(i,j)=1.0;
		
	}
}
}

void solve_x_mom(dmatrix&u,dmatrix&unp1,dmatrix	&v,dmatrix &p)
{

double Ax,Ay,Dx,Dy,gradp_x;

#pragma omp parallel for schedule(static) private(Ax,Ay,Dx,Dy,gradp_x)
for(int i=1;i<=(nx-2);i++)
{
	for(int j=1;j<=(ny-1);j++)
	{
		Ax=(u(i+1,j)*u(i+1,j) - u(i-1,j)*u(i-1,j)) / (2.0*dx);
		
		Ay=0.25 *(((u(i,j) + u(i,j+1)) * (v(i,j) + v(i+1,j)) - (u(i,j) + u(i,j-1)) * (v(i+1,j-1) + v(i,j-1)) ) / dy );
		
		Dx=(u(i+1,j)+u(i-1,j)-2*u(i,j))/(dx*dx);
		
		Dy=(u(i,j+1)+u(i,j-1)-2*u(i,j))/(dy*dy);
		
		gradp_x=p(i+1,j)-p(i,j);
		
		unp1(i,j)=u(i,j)-(Ax+Ay)*dt -(dt/dx)*(gradp_x) +(dt/Re)*(Dx+Dy);
	
	
	}

}

}


void solve_y_mom(dmatrix&v,dmatrix &vnp1,dmatrix &u,dmatrix &p)
{

double Ax,Ay,Dx,Dy,gradp_y;

#pragma omp parallel for schedule(static) private(Ax,Ay,Dx,Dy,gradp_y)
for (int i = 1; i <= nx - 1; i++)
 {
    for (int j = 1; j <= ny - 2; j++)
	 {

		Ax=0.25*(((u(i,j) + u(i,j+1)) * (v(i,j) + v(i+1,j)) - (u(i-1,j) + u(i-1,j+1)) * (v(i,j) + v(i-1,j)))/ dx);

		Ay=(v(i,j+1) * v(i,j+1) - v(i,j-1) * v(i,j-1)) / (2.0 * dy);

		Dx=(v(i+1,j)+v(i-1,j)-2*v(i,j))/(dx*dx);
		
		Dy=(v(i,j+1)+v(i,j-1)-2*v(i,j))/(dy*dy);
		
		gradp_y=p(i,j+1)-p(i,j);
		
		vnp1(i,j)=v(i,j)-(Ax+Ay)*dt -(dt/dy)*(gradp_y) +(dt/Re)*(Dx+Dy);

		
    }
}

}



void solve_pressure(dmatrix&u,dmatrix &v,dmatrix &p,dmatrix&pnp1)
{

#pragma omp parallel for schedule(static)
for(int i = 1; i <= (nx - 1); i++)
{
	for(int j = 1; j <= (ny - 1); j++)
	{
		pnp1(i,j) = p(i,j) - dt * delta * ( (u(i,j) - u(i-1,j)) / dx + (v(i,j) - v(i,j-1)) / dy );
	}
}

}



void apply_bc(dmatrix&u,dmatrix &v,dmatrix &p)
{


//u velocity 

for(int j = 1; j <= ny - 1; j++)
{
	u(0,j) = 0.0;
	u(nx - 1,j) = 0.0;
}

for(int i = 0; i <= nx - 1; i++)
{
	u(i,0) = -u(i,1);
	u(i,ny) = 2 - u(i,ny - 1);
	
}
//v velocity

for(int j = 1; j <= ny - 2; j++)
{
	v(0,j) = -v(1,j);
	v(nx,j) = -v(nx - 1,j);	
}

for(int i = 0; i <= nx; i++)
{
	v(i,0) = 0.0;
	v(i,ny - 1) = 0.0;
}


//pressure
for (int i = 1; i <= (nx - 1); i++)
{
	p(i,0) = p(i,1);
	p(i,ny) = p(i,ny - 1);
}

for (int j = 0; j <= ny; j++)
{
	p(0,j) = p(1,j);
	p(nx,j) = p(nx - 1,j);
}



}


double compute_residual(dmatrix&u,dmatrix &v)
{
double res=0;

double mb;

#pragma omp parallel for schedule(static) private(mb) reduction(+:res)
for (int i = 1; i <= (nx - 1); i++)
{
	for (int j = 1; j <= (ny - 1); j++)
	{
		mb = ( (u(i,j) - u(i-1,j)) / dx + (v(i,j) - v(i,j-1)) / dy );
		res+= fabs(mb);
		
	}
}

return res;

}



void compute_collocated_values(dmatrix &u,dmatrix &v ,dmatrix &p,dmatrix &uc,dmatrix &vc ,dmatrix &pc)
{


for (int i = 0; i <= (nx - 1); i++)
{
	for (int j = 0; j <= (ny - 1); j++)
	{	
		uc(i,j) = 0.5 * (u(i,j) + u(i,j+1));
		vc(i,j) = 0.5 * (v(i,j) + v(i+1,j));
		pc(i,j) = 0.25 * (p(i,j) + p(i+1,j) + p(i,j+1) + p(i+1,j+1));
	}
}

}




void swap_variables(dmatrix &u,dmatrix &v ,dmatrix &p,dmatrix &unp1,dmatrix &vnp1,dmatrix &pnp1)
{

#pragma omp parallel for schedule(static)
for (int i = 0; i <= (nx - 1); i++)
{
	for (int j = 0; j <= ny; j++)
	{
		u(i,j) = unp1(i,j);
	}
}

#pragma omp parallel for schedule(static)
for (int i = 0; i <= nx; i++)
{
	for (int j = 0; j <= (ny - 1); j++)
	{	
		v(i,j) = vnp1(i,j);
	}
}

#pragma omp parallel for schedule(static)

for (int i = 0; i <= nx; i++)
{
	for (int j = 0; j <= ny; j++)
	{
		p(i,j) = pnp1(i,j);
	}
}

}




