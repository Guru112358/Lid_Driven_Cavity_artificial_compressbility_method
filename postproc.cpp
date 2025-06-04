void write_file(dmatrix &u,dmatrix &v ,dmatrix &p ,int nx,int ny);





void write_file(dmatrix &u,dmatrix &v ,dmatrix &p ,int nx,int ny)
{
 	 std::string str="final";
	std::string ext=".csv";

	//str=str.append(std::to_string(number));
	str=str.append(ext);

	std::fstream data;

	data.open(str,std::ios::out);
	data<<" x , y , u , v , p"<<std::endl;
			


for (int i=0; i<nx; i++)
{
	for (int j=0; j<ny; j++)
	
	{
	
        data<<i*dx<<" , "<<j*dy<<" , "<<u(i,j)<<" , "<<v(i,j)<<" , "<<p(i,j)<<std::endl;
					
	}
	
}

    data.close();

}
