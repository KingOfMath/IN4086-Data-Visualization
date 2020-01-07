const config = {
    // Encoding is not recommended to be modified.
    // Instead, it is recommended to copy the contents of the CSV file produced by yourself to example.csv.
    // The encoding format of example.csv is supported by all languages.
    encoding: "UTF-8",

    max_number: 10,

    showMessage: true,

    // Auto Sort by Time
    // Please ensure using standard datetime format (YYYY-MM-DD HH:MM) when this term is enabled!!!
    auto_sort: false,

    timeFormat: "%Y-%m",

    reverse: false,

    divide_by: 'type',

    divide_color_by: 'name',

    color: "blue",

    changeable_color: true,

    divide_changeable_color_by_type: false,
    color_range: ['#830707', '#e13028'],

    // left label
    itemLabel: "Left",

    // right label
    typeLabel: "Right",

    // Top item information horizontal location
    item_x: 50,

    interval_time: 1,

    text_y: -50,

    text_x: 1000,

    offset: 350,

    // Hide barInfo if bar is shorter than barInfo
    display_barInfo: 0,

    use_counter: false,

    step: 1,

    format: ",.0f",

    postfix: "",

    deformat: function (val, postfix) {
        return Number(val.replace(postfix, "").replace(/\,/g, ""));
    },

    left_margin: 0,
    right_margin: 50,
    top_margin: 0,
    bottom_margin: 0,

    dateLabel_switch: true,
    dateLabel_x: null,
    dateLabel_y: null,

    allow_up: false,

    always_up: false,

    enter_from_0: true,

    big_value: false,

    use_semilogarithmic_coordinate: true,

    long: false,

    wait: 0,

    update_rate: 1,

    showLabel: false,

    labelx: -10,

    use_img: false,

    imgs: {
        "条目": "http://i1.hdslb.com/bfs/face/983034448f81f45f05956d0455a86fe0639d6a36.jpg",
        "任意名称": "path/to/img"
    },

    background_color: "#000000",

    rounded_rectangle: true,

    show_x_tick: true,

    // limit bar info display length
    bar_name_max: 30
};