# initialize an empty dictionary to store the aggregate data
aggregates = {}

while True:
    # get input from user
    aggregate_id = input("Enter Aggregate ID (or type 'done' to exit): ")
    
    if aggregate_id == 'done':
        break
    
    dilution = float(input("Enter dilution factor: "))
    media_type = input("Enter media type: ")
    cfu = int(input("Enter CFU: "))
    
    # calculate CFU/g
    cfu_g = (cfu / 0.1) * (1 / dilution)
    
    # add data to dictionary
    if aggregate_id not in aggregates:
        aggregates[aggregate_id] = []
        
    aggregates[aggregate_id].append({
        'dilution': dilution,
        'media_type': media_type,
        'cfu': cfu,
        'cfu_g': cfu_g
    })

# print the table
print("Aggregate ID\tDilution\tMedia Type\tCFU\tCFU/g")
for aggregate_id, data in aggregates.items():
    for datum in data:
        print(f"{aggregate_id}\t{datum['dilution']}\t{datum['media_type']}\t{datum['cfu']}\t{datum['cfu_g']}")

