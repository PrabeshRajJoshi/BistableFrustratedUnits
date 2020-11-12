'''
Python module to execute FORTRAN scripts for simulation of Bistable Frustrated Units (BFUs), to create plots, lattice diagrams and movies, and to organize
the files according to simulation parameters.
'''

def execute_fortran( script_name ,NROWS,BETA_R,T_FINAL,T_STEP,ALPHA_START,ALPHA_STEP
,ALPHA_SPEED,WRITE_STEP,A_SEED,B_SEED,FILENAME) :
    '''
    function to execute the fortran code
    '''
    
    '''
    flist = glob.glob("∗A.dat") ## check for old data files
    if len(flist) > 0:
        print("\nOld data files found!")
        quit("Please remove all data files from this directory!\n")
    '''

    if os.path.exists(script_name+".f90") == False:
        quit("Fortran script does not exist!")
    
    print("Executing Fortran script...\n")
    command_compile_fortran = "gfortran " + script_name+".f90 " + "-o "+script_name+".out"
    option_line1 = " nx "+NROWS + " ny "+NROWS + " bR "+BETA_R + " rt "+T_FINAL + " dt "+T_STEP + " a0 "+ALPHA_START
    option_line2 = " a1 "+ALPHA_SPEED + " sa "+ALPHA_SPEED + " ws "+WRITE_STEP + " sA "+A_SEED + " sB "+B_SEED + " fn "+FILENAME
    execute_options = option_line1 + option_line2
    command_execute_fortran = "time ./" + script_name +".out" + execute_options

    call(command_compile_fortran, shell=True)
    call(command_execute_fortran, shell=True)


def imaging_grid( line, nrow, count, fig_name_add):
    '''
    function to call as multiprocess
    '''
    fig_name = str(count) + fig_name_add

    # make matrix out of the line of data (only for numpy list)
    line = line.reshape((nrow, -1))
    plt.matshow(line)
    plt.savefig(fig_name)


def imaging_and_movie(prefix, write_step, protein_name, movie_decision, stop_time):
    '''
    function for imaging
    '''
    #get the file for imaging
    data_file = prefix + "*" + protein_name + "*.dat"

    # get the name of the file to make image from
    fname = glob.glob(data_file)[0]
    if len(fname)<1:
        print("Data file %s not found!" %data_file)
        quit()
    else:
        print("File used for imaging: ", fname)

    nrows = int(fname[0:3])
    # end part of figure name
    fig_name_end = "_" + prefix + ".png"
    #name of the movie
    movie_name = prefix + ".mp4"

    # starting image number (used to name image files)
    i_num = 1
    # tracks the time depending on "write_step"
    t_count = 0
    print("Executing imaging process...\n")

    ifile = open(fname,'r')
    # iterate over lines in a text file
    for line in ifile:
        if t_count%5 == 0:
            if t_count > stop_time:
                break
            line = np.array( map(float, line.split() ) )
            proc = Process(target=imaging_grid, args=(line, nrows, i_num, fig_name_end))
            proc.daemon = True
            proc.start()
            i_num += 1
        t_count += write_step
    ifile.close()

    if movie_decision == "yes":
        ## make video from images (to install ffmpeg: sudo apt-get install ffmpeg)
        command_make_video = "ffmpeg -qscale 1 -r 20 -b 9600 -i %d" + fig_name_end + " " + movie_name
        call(command_make_video, shell=True)
        print("Finished making movie(frame rate = 20).\n")


def organize_files(dir_name, prefix, script_name, alpha_step, beta_R):
    available = os.path.exists

    print("Organizing files...\n")

    #move data and image/movie file to given directory name "dir_name"
    command = "mkdir "+ dir_name
    call(command, shell=True)

    command = "mv " + prefix + "*.dat" + " " + dir_name
    call(command, shell=True)

    command_move_segment_plot = "mv " + prefix + "*.png" + " " + dir_name
    call(command_move_segment_plot, shell=True)

    if available(prefix + ".mp4"):
        command_remove_images = "rm -rf *_" + prefix + ".png"
        call(command_remove_images, shell=True)

        command_move_mp4 = "mv "+ prefix + ".mp4" + dir_name
        call(command_move_mp4, shell=True)
    
    else:
        command = "mv *_" + prefix + ".png" + dir_name
        call(command, shell=True)
    
    # move "dir_name" directory to "changing_alpha"/"constant_alpha" path
    if float(alpha_step) > 0.0:
        data_path = "changing_alpha/" + script_name + "/bR_" + beta_R
        
    elif float(alpha_step) == 0.0:
        data_path = "constant_alpha/" + script_name + "/bR_" + beta_R
    
    if available(data_path) == False:
        # make required parent directory too
        call("mkdir -p " + data_path, shell=true)
    
    command = "mv " + dir_name + " " + data_path
    call(command, shell=True)
    print("Images/movie and data files moved to directory: \n%s"%data_path)


def segment_plots(prefix, protein_name, n_bfu, write_step, start_time, stop_time, plot_title):
    '''
    function to plot for inherent time
    '''
    A_file = prefix + "*" + protein_name + ".dat"
    # get the name of the file to make the image from
    fname = glob.glob(A_file)[0]
    print("File used for imaging: ",fname)

    n_bfu = int(n_bfu)
    write_step = int(float(write_step))
    start_time = int(start_time)
    stop_time = int(stop_time)
    mid_id = int(n_bfu/2) + 15

    # make the array of evolution of 6 oscillators
    image = []
    t_count = 0
    add_image = image.append
    ifile = open(fname, "r")

    #iterate over lines in text file
    for line in ifile:
        if t_count < start_time:
            t_count += write_step
            continue
        if t_count >= stop_time:
            break

        line = map(float, line.split() )
        ## new shortened array containing only 6 oscillators
        # unit_array = line[1:4] + line[mid_id:mid_id+3]
        # add_image(unit_array)

        add_image(line)
        t_count += write_step

    image = np.array(image)
    image = np.transpose(image)
    ifile.close()


    data_length = len(image[1])
    print("Number of data points for each unit: ", data_length)

    # time array
    x_series = np.arange( start_time, stop_time, save_step )
    print("Number of data points in time array: ", len(x_series))

    x_label = "$time$"
    y_label = "$" + protein_name + "(t)$"

    plt.figure(figsize = (12,7))
    count = 1
    for y_series in image:
        if len(x_series) != len(y_series):
            print("segment_plots: x and y series have different lengths!\n")
            print("x_series: %d y_series: %d", %(len(x_series, len(y_series))))
            quit()
        else:
            plt.plot(x_series, y_series, label = str(count))
            count += 1
    # add in labels and title
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(plot_title)
    plt.grid(True)
    # plt.legend()

    # save image as PNG file
    plt.savefig(prefix + "_" + protein_name + ".png", dpi=300)
    plt.close()