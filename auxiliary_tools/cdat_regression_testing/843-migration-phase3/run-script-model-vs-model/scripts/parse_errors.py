import re


def parse_unique_errors_without_timestamp(log_file_path):
    unique_errors = set()
    timestamp_pattern = re.compile(r"^\d{2}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} ")

    with open(log_file_path, "r") as file:
        for line in file:
            if "Error" in line or "error" in line or "ERROR" in line:
                # Remove the timestamp
                cleaned_line = timestamp_pattern.sub("", line).strip()
                unique_errors.add(cleaned_line)

    print("Unique error messages without timestamp:")
    for error in unique_errors:
        print(error)


# Replace '24-09-24-main-log.txt' with the path to your log file
# log_file_path = "auxiliary_tools/cdat_regression_testing/843-migration-phase3/run-script-model-vs-model/24-09-24-main-log.txt"
log_file_path = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/843-migration-phase3/run-script-model-vs-obs/25-09-24-log.txt"
parse_unique_errors_without_timestamp(log_file_path)
