/**
 * dump calculated bin files for quick checking
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-20
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <fstream>
#include <iostream>


int main(int argc, char** argv)
{
	if(argc<=1)
	{
		std::cerr << "Give a bin file" << std::endl;
		return -1;
	}

	std::ifstream ifstr(argv[1], std::ios_base::binary);
	if(!ifstr)
	{
		std::cerr << "Cannot read " << argv[1] << std::endl;
		return -1;
	}


	while(!ifstr.eof())
	{
		unsigned int num = 0;
		ifstr.read((char*)&num, sizeof(num));
		std::cout << "Number of branches: " << num << std::endl;

		for(unsigned int i=0; i<num; ++i)
		{
			double E=0, w1=0, w2=0, w3=0;
			ifstr.read((char*)&E, sizeof(E));
			ifstr.read((char*)&w1, sizeof(w1));
			ifstr.read((char*)&w2, sizeof(w2));
			ifstr.read((char*)&w3, sizeof(w3));

			std::cout << "#" << i << ": " << "E=" << E << "ws=[" << w1 << ", " << w2 << ", " << w3 << "]" << std::endl;
		}
	}


	return 0;
}
