using CSV

#Choices for Description : Snf1:glucose_repression | wt:cAMP | sch9Delta:cAMP | Sch9P:100_glutamine | Sch9P:1_glutamine |
# Sch9:gtr1Delta | Sch9:glucose_starve | Sch9:glucose_relief | Mig1:glucose_relief | Rib:rap | Gln3:rap

#This code gives a dictionary with time points, experimentall values and pre shift and post shift values of
#some parameters like ex: Carbon and ATP.

#Exampel of use:
#Description = Sch9:glucose_starve
#Dictionary = get_data_from_CSV(Description)
#Time_points = Dictionary["Time_points"]
#Experimental_values = Dictionary["Experimental_values"]
#Pre_shift_carbon = Dictionary["pre_shift"]["Carbon"]

function get_data_from_CSV(Description)

    function looking_for_row_number(Description)
        row_number = 0
        for row in p["Description"]
            row_number += 1
            #print(row)
            if Description == row
                break
            end
        end
        return row_number
    end

    function get_time_points(row_number, Description)
        #Creates a substring
        time_list = split(p["Time points"][row_number],",")
        N = 0
        #Removes some characters from each substring
        for item in time_list
            N += 1
            removechar = ['[', ']']
            time_list[N] = replace(item, removechar => "")
        end
        #Creates list of all numbers that now can be used (string => float)
        Time_points = parse.(Float64,time_list)
        # These desriptions have s instead of min, so I divide by 60 on these to get the same unit.
        if Description == "Mig1:glucose_relief" || "sch9Delta:cAMP" || "wt:cAMP"
            Time_points = Time_points/60
        end
        return Time_points
    end

    function get_experimental_points(row_number)
        #Creates a substring
        exp_list = split(p["Experimental values"][row_number],",")
        N = 0
        #Removes some characters from each substring
        for item in exp_list
            N += 1
            removechar = ['[', ']']
            exp_list[N] = replace(item, removechar => "")
        end
        #Creates list of all numbers that now can be used (string => float)
        Experimental_values = parse.(Float64,exp_list)
        return Experimental_values
    end

    function looking_for_pre_shift_variable_values_and_names(row_number)
        #Creates a substring
        pre_shift_list = split(p["Simulation pre shift specification"][row_number],(',',':',''','}','{'))
        pre_shift_list = [i for i in pre_shift_list if i != ""]
        pre_shift_list = [i for i in pre_shift_list if i != "inconds"]
        pre_shift_list = [i for i in pre_shift_list if i != "parameters"]
        pre_shift_values = []
        pre_shift_parameters = pre_shift_list
        for item in pre_shift_list
            try
                value = parse.(Float64,item)
                append!(pre_shift_values,value)
                pre_shift_parameters = [i for i in pre_shift_parameters if i != item]
            catch
                continue
            end
        end
        return pre_shift_parameters, pre_shift_values
    end

    function looking_for_post_shift_variable_values_and_names(row_number)
        #Creates a substring
        post_shift_list = split(p["Simulation post shift specification"][row_number],(',',':',''','}','{'))
        post_shift_list = [i for i in post_shift_list if i != ""]
        post_shift_list = [i for i in post_shift_list if i != "inconds"]
        post_shift_list = [i for i in post_shift_list if i != "parameters"]
        post_shift_values = []
        post_shift_parameters = post_shift_list
        for item in post_shift_list
            try
                value = parse.(Float64,item)
                append!(post_shift_values,value)
                post_shift_parameters = [i for i in post_shift_parameters if i != item]
            catch
                continue
            end
        end
        return post_shift_parameters, post_shift_values
    end

    function create_pre_shift_dictionary(pre_shift_parameters, pre_shift_values)
        pre_shift_dictionary = Dict{String, Float64}()
        N = 0
        for i in pre_shift_parameters
            N += 1
            merge!(pre_shift_dictionary,Dict(pre_shift_parameters[N] => pre_shift_values[N]))
        end
        return pre_shift_dictionary
    end

    function create_post_shift_dictionary(post_shift_parameters, post_shift_values)
        post_shift_dictionary = Dict{String, Float64}()
        N = 0
        for i in post_shift_parameters
            N += 1
            merge!(post_shift_dictionary,Dict(post_shift_parameters[N] => post_shift_values[N]))
        end
        return post_shift_dictionary
    end

    ######MAIN#######

    p = CSV.File(pwd()*"/Data/mc-e20-02-0117-s04.CSV")

    row_number = looking_for_row_number(Description)

    Time_points = get_time_points(row_number, Description)

    Experimental_values = get_experimental_points(row_number)

    pre_shift_parameters, pre_shift_values = looking_for_pre_shift_variable_values_and_names(row_number)

    post_shift_parameters, post_shift_values = looking_for_post_shift_variable_values_and_names(row_number)

    pre_shift_dictionary = create_pre_shift_dictionary(pre_shift_parameters, pre_shift_values)

    post_shift_dictionary = create_post_shift_dictionary(post_shift_parameters, post_shift_values)

    Dictionary_of_parameters_and_values = Dict("Time_points" => Time_points, "Experimental_values" => Experimental_values, "pre_shift" => pre_shift_dictionary, "post_shift" => post_shift_dictionary)
    
    return Dictionary_of_parameters_and_values
end





