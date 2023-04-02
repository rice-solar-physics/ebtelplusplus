// ****
// *
// * Routine that uses TinyXML library to parse XML configuration files
// *
// * Author: Will Barnes
// *
// * Date last modified: 16/12/2015
// *
// ****

#include "xmlreader.h"

std::string get_element_text(tinyxml2::XMLElement *root,std::string search_value)
{
	//Search for element in root
	tinyxml2::XMLElement * result = recursive_read(root,search_value);
	tinyxml2::XMLElement * checked_result = check_element(result,search_value);
	return checked_result->GetText();
}

tinyxml2::XMLElement * get_element(tinyxml2::XMLElement *root,std::string search_value)
{
	tinyxml2::XMLElement * result = recursive_read(root,search_value);
	tinyxml2::XMLElement * checked_result = check_element(result,search_value);
	return checked_result;
}

tinyxml2::XMLElement * check_element(tinyxml2::XMLElement *result,std::string search_value)
{
	if(result != NULL)
	{
		return result;
	}
	else
	{
		//TODO:default value when value cannot be found; search for default configuration.
		//Package default XML with HYDRAD
		printf("Could not find %s.\n",search_value.c_str());
		printf("Returning NULL.\n");
		return NULL;
	}
}

tinyxml2::XMLElement * recursive_read(tinyxml2::XMLElement *parent, std::string search_value)
{
	tinyxml2::XMLElement *child;

	if(search_value.compare(parent->Value()) == 0)
	{
		child = parent;
	}
	else
	{
		child = parent->FirstChildElement();
	}

	while( child != NULL && search_value.compare(child->Value()) != 0)
	{
		if(child->FirstChildElement() != NULL)
		{
			tinyxml2::XMLElement *grandChild = recursive_read(child,search_value);
			if(grandChild != NULL && search_value.compare(grandChild->Value()) == 0)
			{
				child = grandChild;
				break;
			}
		}

		child = child->NextSiblingElement();
	}

	return child;
}

bool string2bool(std::string boolString)
{
	bool boolReturn;
	
	if(boolString.compare("true")==0 || boolString.compare("True")==0)
	{
		boolReturn = true;
	}
	else if(boolString.compare("false")==0 || boolString.compare("False")==0)
	{
		boolReturn = false;
	}
	else
	{
		printf("Unrecognized boolean option. Defaulting to false.\n");
		boolReturn = false;
	}
	
	return boolReturn;
}