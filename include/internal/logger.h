/*
 * logger.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef LOGGER_H_
#define LOGGER_H_

#include <string>
#include <sstream>
#include <iostream>
#include "termcolor.hpp"

namespace FastVQA{

#define LOG_LEVEL 1 //0 - debug, 1 - info, 2 - warning, 3 - error

inline std::string getCurrentDateTime(){

    time_t now = time(0);
    struct tm  tstruct;
    char  buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return std::string(buf);

};
template<typename Functor>
inline void printLog(std::string logType, Functor termColor, std::string logMsg){

	std::string now = getCurrentDateTime();
	std::cerr << termColor << "[[" << now << "]][" << logType << "] " << logMsg << termcolor::reset << std::endl;

}

inline void logd( std::string logMsg, int log_level = LOG_LEVEL){

	if(log_level == 0)
		printLog("DEBUG", termcolor::cyan, logMsg);

}

inline void logi( std::string logMsg, int log_level = LOG_LEVEL){

	if(log_level <= 1)
		printLog("INFO", termcolor::white, logMsg);

}

inline void logw( std::string logMsg, int log_level = LOG_LEVEL){

	if(log_level <= 2)
		printLog("WARNING", termcolor::yellow, logMsg);

}

inline void loge( std::string logMsg ){

	printLog("ERROR", termcolor::red, logMsg);

}

inline void throw_runtime_error( std::string logMsg ){

	printLog("ERROR", termcolor::red, logMsg);
	throw std::runtime_error(logMsg);

}

}

#endif /* LOGGER_H_ */
