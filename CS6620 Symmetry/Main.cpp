#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "pmvs_base/pmvs/patch.h"

void main()
{
	std::vector<Patch::Cpatch> patches;

	for(int i = 0; i < 9; i++)
	{
		// file names expected to be 0000.patc0, 0001.patc0, etc.
		std::stringstream filename;
		filename << std::setw(4) << std::setfill('0') << i << ".patc0";

		std::cout << "Starting to read image " << filename.str() << std::endl;

		std::ifstream data;
		data.open(filename.str().c_str());

		if (data.is_open())
		{
			std::string line;

			// value of "A" from PMVS-1's match.exe (ignored)
			getline(data, line);

			// number of points in file
			int points;
			data >> points;
			
			std::cout << "Points in image " << i << ": " << points << std::endl;

			//int count = 0;
			while (data.good())//&& count < points)
			{
				// skip lines until find "PATCH0"
				//std::streampos start;
				do
				{
					//start = data.tellg();
					getline(data, line);

					//std::cout << "line is " << line << std::endl;
				}
				while (data.good() && !(line == "PATCH0"));

				if (data.good())
				{
					//data.seekg(start); // for some reason this doesn't work

					// read in
					Patch::Cpatch p;
					data >> p;
					patches.push_back(p);
					
					//std::cout << p << std::endl;
				}
				//count++;
			}
			data.close();
		}
		else
		{
			std::cout << "Couldn't read file!" << std::endl;
		}
	}
	std::cout << "Done! enter some number to exit..." << std::endl;
	int i;
	std::cin >> i;
}
