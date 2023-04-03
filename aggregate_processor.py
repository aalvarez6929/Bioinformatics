import csv

# define ANSI escape codes for text colors
color_map = {
    'GYE': '\033[38;5;202m',   # orange
    'R2A': '\033[38;5;99m',    # purple
    'Cellulose': '\033[38;5;40m',  # green
    'VL55': '\033[38;5;27m'    # blue
}

# initialize an empty dictionary to store the aggregate data and photo names
aggregates = {}
photo_names = {}

# function to calculate dilution
def calculate_dilution(mass, volume):
    return mass / volume

# function to add data to dictionary
def add_data(aggregate_id, mass, volume, media_type, cfu, photo_id):
    # calculate dilution
    dilution = 0
    if volume > 0:
        dilution = calculate_dilution(mass, volume)
    
    # calculate CFU/g
    cfu_g = 0
    if dilution > 0:
        cfu_g = (cfu / 0.1) * (1 / dilution)
    
    # add data to dictionaries
    if aggregate_id not in aggregates:
        aggregates[aggregate_id] = []
        
    aggregates[aggregate_id].append({
        'mass': mass,
        'volume': volume,
        'media_type': media_type,
        'cfu': cfu,
        'dilution': dilution,
        'cfu_g': cfu_g,
        'photo_id': photo_id
    })
    
    if aggregate_id not in photo_names:
        photo_names[aggregate_id] = []
    
    photo_names[aggregate_id].append(photo_id)

# read data from input file
input_filename = input("Enter input filename: ")
with open(input_filename, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    headers = next(reader)
    for row in reader:
        try:
            aggregate_id, mass, volume, media_type, cfu, photo_id = row
            add_data(aggregate_id, float(mass), float(volume), media_type, int(cfu), photo_id)
        except ValueError:
            # handle missing values
            pass

# read photo names from separate file
photo_filename = input("Enter photo filename: ")
with open(photo_filename, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    next(reader)  # skip header row
    for row in reader:
        try:
            aggregate_id, photo_names_str = row
            photo_names[aggregate_id] = photo_names_str.split(', ')
        except ValueError:
            # handle missing values
            pass

# write data to output file
output_filename = input("Enter output filename: ")
with open(output_filename, 'w') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(['Aggregate ID', 'Mass (g)', 'Volume (mL)', 'Media Type', 'Photo ID', 'CFU', 'Dilution', 'CFU/g'])
    for aggregate_id, data in aggregates.items():
        for datum in data:
            media_type = datum['media_type']
            writer.writerow([aggregate_id, datum['mass'], datum['volume'], media_type, datum['photo_id'], datum['cfu'], '{:.2e}'.format(datum['dilution']), '{:.2e}'.format(datum['cfu_g'])])

## write photo names to separate file
photo_filename = input("Enter photo filename: ")
with open(photo_filename, 'w') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(['Aggregate ID', 'Photo Names'])
    for aggregate_id, photo_names_list in photo_names.items():
        photo_names_str = ', '.join(photo_names_list)
        writer.writerow([aggregate_id, photo_names_str])
