/*
PowerLaw Radiation Tests
Different tests for radiation model as applied here
*/

#include "boost/program_options.hpp"
#include "../source/helper.h"
#include "../Radiation_Model/source/radiation.h"

int main(int argc, char *argv[])
{
  // Command line arguments
  namespace po = boost::program_options;
  po::options_description description("Tests for ebtel++");
  description.add_options()
    ("help,h","The help message")
    ("quiet,q",po::bool_switch()->default_value(false),"Suppress output.")
    ("output_file,o",po::value<std::string>()->default_value("test.out"),"Output file from tests.");
  po::variables_map vm;
  po::store(po::command_line_parser(argc,argv).options(description).run(), vm);
  if(vm.count("help"))
  {
  	std::cout << description;
  	return 0;
  }
  po::notify(vm);

  // Declare Radiation model
  PRADIATION radiation_model;
  radiation_model = new CRadiation();
  // Make temperature vector over relevant range
  double temperature_min = pow(10.0,4.0);
  double temperature_max = pow(10.0,8.0);
  std::vector<double> temperature (1000);
  double delta_temperature = (log10(temperature_max) - log10(temperature_min))/temperature.size();
  std::vector<double> radiation (1000);

  // Call radiation model
  for(int i=0;i<temperature.size();i++)
  {
    temperature[i] = temperature_min*pow(10.0,i*delta_temperature);
    radiation[i] = radiation_model->GetPowerLawRad(log10(temperature[i]));
  }

  // Print output
  std::ofstream f;
  f.open(vm["output_file"].as<std::string>());
  for(int i=0;i<temperature.size();i++)
  {
    f << temperature[i] << "\t" << radiation[i] << "\n";
  }
  f.close();

  delete radiation_model;

  return 0;
}
