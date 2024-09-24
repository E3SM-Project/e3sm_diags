def parse_unique_errors(log_file_path):
    unique_errors = set()

    with open(log_file_path, "r") as file:
        for line in file:
            if "OSError" in line:
                unique_errors.add(line.strip())

    print("Unique error messages:")
    for error in unique_errors:
        print(error)


# Replace '24-09-24-main.log' with the path to your log file
log_file_path = "auxiliary_tools/cdat_regression_testing/843-migration-phase3/run-script-model-vs-model/24-09-24-main.log"
parse_unique_errors(log_file_path)
