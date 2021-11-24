import sys
import datetime

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
YELLOW = '\033[93m'
CYAN = "\033[36m"
BLUE = "\033[34m"
DIM = '\033[2m'


def colour(text, text_colour):
    bold_text = 'bold' in text_colour
    text_colour = text_colour.replace('bold', '')
    underline_text = 'underline' in text_colour
    text_colour = text_colour.replace('underline', '')
    text_colour = text_colour.replace('_', '')
    text_colour = text_colour.replace(' ', '')
    text_colour = text_colour.lower()
    if 'red' in text_colour:
        coloured_text = RED
    elif 'green' in text_colour:
        coloured_text = GREEN
    elif 'yellow' in text_colour:
        coloured_text = YELLOW
    elif 'blue' in text_colour:
        coloured_text = BLUE
    elif 'cyan' in text_colour:
        coloured_text = CYAN
    elif 'magenta' in text_colour:
        coloured_text = MAGENTA
    elif 'dim' in text_colour:
        coloured_text = DIM
    else:
        coloured_text = ''
    if bold_text:
        coloured_text += BOLD
    if underline_text:
        coloured_text += UNDERLINE
    if not coloured_text:
        return text
    coloured_text += text + END_FORMATTING
    return coloured_text



def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def print_log(text):
    sys.stderr.write(colour(get_timestamp(),"bold")+" "+colour(text,"bold_yellow")+"\n")
    sys.stderr.flush()

def print_message(text):
    sys.stderr.write("   "+colour(text,"green")+"\n")
    sys.stderr.flush()

def print_error(text, error_number = 1, exit = True):
    sys.stderr.write("   "+colour("Error: ","red_bold")+colour(text,"red")+"\n")
    sys.stderr.flush()
    if exit:
        if error_number != 0:
            sys.exit(error_number)
        else:
            sys.exit(1)

def print_warning(text):
    sys.stderr.write("   "+colour("Warning: ","magenta_bold")+colour(text,"magenta")+"\n")
    sys.stderr.flush()
