/*
 * logger.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LOGGER_H_
#define SRC_LOGGER_H_

#include <string>
#include <sstream>
#include <iostream>

#define LOG_LEVEL 0 //0 - debug, 1 - info, 2 - warning, 3 - error

inline std::string getCurrentDateTime(){

    time_t now = time(0);
    struct tm  tstruct;
    char  buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return std::string(buf);

};

inline void printLog(std::string logType, std::string logMsg){

	std::string now = getCurrentDateTime();
	std::cout << "[[" << now << "]][" << logType << "] " << logMsg << std::endl;

}

inline void logd( std::string logMsg ){

	if(LOG_LEVEL == 0)
		printLog("DEBUG", logMsg);

}

inline void logi( std::string logMsg ){

	if(LOG_LEVEL <= 1)
		printLog("INFO", logMsg);

}

inline void logw( std::string logMsg ){

	if(LOG_LEVEL <= 2)
		printLog("WARNING", logMsg);

}

inline void loge( std::string logMsg ){

	printLog("ERROR", logMsg);

}


#endif /* SRC_LOGGER_H_ */
