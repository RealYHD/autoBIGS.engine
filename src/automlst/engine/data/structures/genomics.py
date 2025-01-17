from dataclasses import dataclass
from numbers import Number
from typing import Mapping, Sequence, Set, Union


@dataclass(frozen=True)
class StringAnnotation:
    type: str
    start: int
    end: int
    feature_properties: Mapping[str, Set[str]]

@dataclass(frozen=True)
class NamedString:
    name: str
    sequence: str

@dataclass(frozen=True)
class AnnotatedString(NamedString):
    annotations: Sequence[StringAnnotation]

@dataclass(frozen=True)
class SangerTraceData(NamedString):
    seq_param_file_name: str
    analysis_proto_settings_name: str
    analysis_rpto_settings_ver: str
    analysis_proto_xml_data: str
    analysis_proto_xml_schema_ver: str
    sample_comment: Union[None, str]
    capillary_machine: bool
    container_identifier: str
    container_name: str
    comment_title: str
    channel_1: Sequence[Number]
    channel_2: Sequence[Number]
    channel_3: Sequence[Number]
    channel_4: Sequence[Number]
    measured_voltage_dv: Sequence[Number]
    measured_current_ma: Sequence[Number]
    measured_power_mw: Sequence[Number]
    measured_temperature_celsius: Sequence[Number]
    down_sample_factor: Number
    dye_1: str
    dye_2: str
    dye_3: str
    dye_4: str
    dye_wavelength_1: str
    dye_wavelength_2: str
    dye_wavelength_3: str
    dye_wavelength_4: str
    dye_set_name: str
    electrophoresis_voltage_setting_v: Number
    start_run_event: str
    stop_run_event: str
    start_collection_event: str
    stop_collection_event: str
    base_order: Sequence[str]
    gel_type_desc: str
    injection_time_sec: Number
    inection_voltage_v: Number
    lane_or_capillary: Number
    sample_tracking_id: str
    length_to_detector_cm: Number
    laser_power_mw: Number
    instrument_name_and_serial: str
    data_collection_module_file: str
    model_number: str
    pixels_avg_per_lane: Number
    number_of_capillaries: Number
    marked_off_scale_scans: Union[None, Sequence[Number]]
    # Skipped Ovrl, OvrV
    mobility_file: str
    # Skipped PRJT, PROJ
    pixel_bin_size: Number
    # Skipped scan rate
    results_group_comment: Union[None, str]
    results_group_name: str
    run_module_ver: str
    run_module_xml: str
    run_module_xml_ver: str
    run_proto_name: str
    run_proto_ver: str
    run_start_date: str  # Date time object
    run_stop_date: str  # Date time object
    data_collection_start_date: str
    data_collection_stop_date: str
    run_name: str
    run_start_time: str  # time object
    run_stop_time: str  # time object
    collection_start_time: str  # time object
    collection_stop_time: str  # time object
    saturated_data_points: Union[None, Sequence[Number]]
    color_rescaling_divisor: Number
    scan_count: Number
    polymer_lot_expiration: str  # date time object
    polymer_lot_number: Number
    sample_name: str
    # Skipped genescan data
    # Skipped size standard file name
    data_collection_software_ver: str
    data_collection_firmware_ver: str
    run_temperature_setting_celcius: Number
    well_id: str
    plate_user_name: str
