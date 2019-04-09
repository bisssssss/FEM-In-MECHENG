cmake_minimum_required (VERSION 2.8.2)
project (main)
add_executable(main main.cc)
add_custom_target(run
	COMMAND main
	DEPENDS main
	${CMAKE_PROJECT_DIR}
	)
