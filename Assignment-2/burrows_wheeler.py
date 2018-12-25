

def burrows_wheeler_transform(string):
    # assert string[-1] == "$", "String must end with `$` to compute BW-Transform"
    # Compute rotations
    len_string = len(string)
    strings = [string[i:] + string[:i] for i in range(len_string, 0, -1)]

    # Sort strings
    strings.sort()

    # Take last column
    bw_transform = [s[-1] for s in strings]

    return "".join(bw_transform)


def inverse_burrows_wheeler_transform(bw_string):

    # Make char list to construct original string
    char_list = [(char, i) for i, char in enumerate(bw_string)]

    # Need to use sorted character list
    sorted_char_list = sorted(char_list)

    # Next char dictionary
    next_char_dict = {char: next_char for char, next_char in zip(sorted_char_list, char_list)}

    next_char = sorted_char_list[0] # It should be `$`

    original = []
    for _ in range(len(char_list)):
        original.append(next_char)
        next_char = next_char_dict[next_char]

    original_str = "".join([c[0] for c in original[::-1]])

    return original_str


def read_transform_type():

    msg = """Please enter the transform type:
    1. BW-Transform 
    2. Inverse BW-Transform
(Enter 1 or 2, Please enter q to quit)
"""

    inp = input(msg)

    # assert isinstance(inp, int), "Input MUST be integer"
    # assert inp == "1" or inp == "2", "You MUST enter 1 or 2"

    if inp == "1":
        return "BW"
    elif inp == "2":
        return "IBW"
    else:
        return inp


def read_input_BW():
    msg = """Enter a string ending with `$` to compute BW-transform:
"""
    inp = input(msg)

    while inp[-1] != "$":
        msg = """Enter a string ending with `$` to compute BW-transform:
"""
        inp = input(msg)

    return inp


def read_input_IBW():
    msg = """Enter a string to compute inverse BW-transform:
"""
    inp = input(msg)

    return inp

if __name__ == "__main__":

    while True:
        transform_type = read_transform_type()
        if transform_type == "BW":
            string = read_input_BW()
            print(burrows_wheeler_transform(string))
        elif transform_type == "IBW":
            string = read_input_IBW()
            print(inverse_burrows_wheeler_transform(string))
        elif transform_type == "q":
            break  
