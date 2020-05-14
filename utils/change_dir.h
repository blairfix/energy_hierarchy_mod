//File: change_dir.h
#ifndef CHANGE_DIR_H
#define CHANGE_DIR_H

#include <stdio.h>
#include <unistd.h>
#include <string.h>

// changes the directory 

void change_directory(std::string rmv, std::string add){

    char*  d = get_current_dir_name();

    //std::cout << d << std::endl;

    std::string dir(d);
    dir.erase(dir.find(rmv));
    dir.append(add);
    const char* dir_new = dir.c_str();

    int ret = chdir(dir_new);
	if (ret != 0){
		perror("chdir() failed");
	}	
	
	
}


#endif
