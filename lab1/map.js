// map GeoJsonLayer rendering
// map bar chart animation
// SPLOM graph

/**
 * draw SPLOM with top 7 attributes
 * @constructor
 */
let constOptions = 0;
function SPLOM() {

    $("svg").empty();

    var rows = [];
    var city = $("#city option:selected").text();

    // filter data
    var mainCityArray = csv.filter(function (v, i) {
        return ((v["name"] === city));
    });

    // make options constant avoiding reloading
    if (constOptions === 0) {
        constOptions = options;
    } else {
        options = constOptions;
    }

    // filter the data with appointed attributes
    options.forEach(function (option) {
        var optionArray = mainCityArray.filter(function (v, i) {
            return ((v["type"] === option));
        });
        let crimeArray = optionArray.map(data => data.crimerate);
        rows.push(crimeArray);
    });

    // add number and set it as value
    var rowsValue = [];
    for (let i = 0; i < options.length; i += 1) {
        const num = rows[i].reduce((a, b) => Number(a) + Number(b), 0);
        rowsValue.push({
            'key': options[i],
            'value': num,
            'row': rows[i]
        });
    }
    rowsValue.sort(function (a, b) {
        return b.value - a.value
    });

    // get key and value
    options = rowsValue.map(data => data.key);
    rows = rowsValue.map(data => data.row);
    let attributes = ["1", "2", "3", "4", "5", "6", "7"];
    options = options.slice(1, attributes.length + 1);
    rows = rows.slice(1, attributes.length + 1);
    colors = d3old.scale.ordinal().range(["#827abf", "#f62150", "#6f89b6", "#f5e0b7", "#5b1e37", "#b9e3c5"]);

    // add a crime type legend
    $("#legend_corr").html("");
    document.getElementById("legend_corr_inner").style.display = 'inline';
    for (var i = 0; i < options.length; i++) {
        var opt = options[i];
        var node = document.createElement("LI");
        var textnode = document.createTextNode(attributes[i] + " = " + options[i]);
        node.appendChild(textnode);
        document.getElementById("legend_corr").appendChild(node);
    }

    //create an n-by-n matrix based on pairs of attributes
    let attributeMatrix = [];
    attributes.forEach(function (a, x) {
        attributes.forEach(function (b, y) {
            attributeMatrix.push({a: a, b: b, x: x, y: y})
        })
    });

    // process the data to be scattered
    var scatterData = new Array(rows[0].length);
    for (let i = 0; i < scatterData.length; i += 1) {
        scatterData[i] = {}
    }
    for (let i = 0; i < rows.length; i += 1) {
        for (let j = 0; j < rows[0].length; j += 1) {
            scatterData[j][attributes[i]] = Number(rows[i][j])
        }
    }

    // read Province data and draw SPLOM
    // create scales dynamically for each attribute's extent
    let scale = {};
    attributes.forEach(function (att) {
        scale[att] = d3old.scale.linear();
        let attExtent = d3old.extent(scatterData, function (d) {
            return d[att];
        });
        scale[att].domain(attExtent).range([5, 75]);
    });

    //bind the matrix array to a grid of g elements
    d3old.select("svg")
        .selectAll("g")
        .data(attributeMatrix)
        .enter()
        .append("g")
        .attr("transform", function (d) {
            return "translate(" + (d.x * 80) + "," + (d.y * 80) + ")"
        });
    d3old.selectAll("g")
        .each(function (matrix, i) {

            //index i is only used for coloring
            //background/border
            d3old.select(this).append("rect").style("fill", "white").style("stroke", "black").style("stroke-width", 1)
                .attr("height", 80)
                .attr("width", 80);

            //label
            d3old.select(this).append("text")
                .attr("x", 40)
                .style("text-anchor", "middle")
                .attr("y", 12)
                .text(matrix.a + " - " + matrix.b);

            //scatter points
            d3old.select(this).selectAll("circle")
                .data(scatterData)
                .enter()
                .append("circle")
                .attr("r", 2)
                .style("fill", colors(i))
                .attr("cx", function (d) {
                    return scale[matrix.a](d[matrix.a])
                })
                .attr("cy", function (d) {
                    return 75 - scale[matrix.b](d[matrix.b])
                })
        });

    // set width and height of the SPLOM
    var svg = d3old.select("#svgrank").append("svg")
        .attr("width", 600)
        .attr("height", 560);
}


/**
 * format the default date type to 'month year'
 * @param date
 * @returns {string}
 */
function formatDate(date) {
    var monthNames = [
        "January", "February", "March",
        "April", "May", "June", "July",
        "August", "September", "October",
        "November", "December"
    ];
    var monthIndex = date.getMonth();
    var year = date.getFullYear();
    return monthNames[monthIndex] + ' ' + year;
}


/**
 * set configs of the time slider
 */
$(function () {
    $("#slider-range").slider({
        range: true,

        // set the time range
        min: new Date('2012-01-01T00:00:00').getTime(),
        max: new Date('2019-11-01T00:00:00').getTime(),
        step: 86400000,

        // set the default range in time slider
        values: [new Date('2012-12-01T00:00:00').getTime(), new Date('2019-11-01T00:00:00').getTime()],

        // set the slider date
        slide: function (event, ui) {
            $("#amount").val(formatDate(new Date(ui.values[0])) + '-' + formatDate(new Date(ui.values[1])));
        },

        // when the animation stops, update data again
        stop: function () {
            filterdata();
        }
    });

    // show date under the time slider
    $("#amount").val(formatDate((new Date($("#slider-range").slider("values", 0)))) +
        " - " + formatDate((new Date($("#slider-range").slider("values", 1)))));
});


/**
 * animate the change of crime rate on the map
 * it is executed recursively
 * @param steps
 * @param currentDate
 */
function AnimationLoop(steps, currentDate) {

    // If the calculated step is zero, then stop animation
    if (steps > 0) setTimeout(function () {

        // every time the date moves on with one month
        currentDate.setMonth(currentDate.getMonth() + 1);
        $("#slider-range").slider('values', 0, currentDate);
        sliderdate = currentDate;

        // update the date under the time slider
        $("#amount").val(formatDate((new Date($("#slider-range").slider("values", 0)))) +
            " - " + formatDate((new Date($("#slider-range").slider("values", 1)))));

        // update the end if stopping date changes, then we can control when to stop
        var stop = new Date($("#slider-range").slider("values", 1));
        steps = (stop.getFullYear() - currentDate.getFullYear()) * 12 + stop.getMonth() - currentDate.getMonth();
        filterdata();

        // recursion
        AnimationLoop(steps, currentDate);
    }, 1000);

};


/**
 * set configs of the animation
 * @constructor
 */
function Animate() {
    var start = new Date($("#slider-range").slider("values", 0));
    var stop = new Date($("#slider-range").slider("values", 1));
    var steps = (stop.getFullYear() - start.getFullYear()) * 12 + stop.getMonth() - start.getMonth();
    var currentDate = start;
    AnimationLoop(steps, currentDate);
}

// define variables
var csv;
var data;
var options;
var sliderdate;
var slidercategory;

// set the configs of the map
Deck = new deck.DeckGL({
    latitude: 52.0907,
    longitude: 5.1214,
    zoom: 7,
    maxZoom: 16,
    pitch: 45,
    container: 'deckgl-container',
    mapboxApiAccessToken: 'pk.eyJ1IjoiZW5qYWxvdCIsImEiOiJjaWhtdmxhNTIwb25zdHBsejk0NGdhODJhIn0.2-F2hS_oTZenAWc0BMf_uw',
    mapStyle: 'mapbox://styles/mapbox/dark-v9'
    //  layers: [geojsonLayer]

});


/**
 * make sure the dataset and json file are completely loaded
 * then function the selection of different crime types
 */
Promise.all([
    d3.json('https://www.webuildinternet.com/articles/2015-07-19-geojson-data-of-the-netherlands/townships.geojson'),
    dataSet
]).then(([townships, crimedata]) => {

    console.log('data loaded!');
    data = townships.features;
    csv = crimedata.map(function (d) {
        return {
            "crimerate": d.value, "name": d.name,
            "period": d.date, "type": d.type
        };
    });
    var select = document.getElementById("dropmenu");
    var temparray = csv.filter(function (v, i) {
        return ((v["period"] === '2012-02' || v["period"] === '2019-02'));
    });

    // filter all the types and return the unique type
    let result = temparray.map(data => data.type);
    options = result.filter(onlyUnique);

    // apply the crime types to the selection dropdown menu
    for (var i = 0; i < options.length; i++) {
        var opt = options[i];
        var el = document.createElement("option");
        el.textContent = opt;
        el.value = opt;
        select.appendChild(el);
    }
    filterdata();
});


/**
 * return the unique type
 * @param value
 * @param index
 * @param self
 * @returns {boolean}
 */
function onlyUnique(value, index, self) {
    return self.indexOf(value) === index;
}


/**
 * select suitable legend due to different scales of crime rate
 */
function legend_select() {
    if ($("#dropmenu :selected").text() === 'Totaal misdrijven') {
        document.getElementById("large").style.display = 'inline';
        document.getElementById("small").style.display = 'none';
    } else {
        document.getElementById("large").style.display = 'none';
        document.getElementById("small").style.display = 'inline';
    }
}


/**
 * update the data on the map according to the date on the time slider
 */
function filterdata() {

    legend_select();

    // get the information of the crime type and date
    sliderdate = new Date($("#slider-range").slider("values", 0));
    slidercategory = $("#dropmenu :selected").text();
    console.log(slidercategory);
    var year = sliderdate.getFullYear();
    var month = sliderdate.getMonth() + 1;
    if (month < 10) {
        month = '0' + month;
    }
    var yearmonth = year + '-' + month;
    console.log(yearmonth);

    var temparray = csv.filter(function (v, i) {
        return ((v["period"] === yearmonth && v["type"] === slidercategory));
    });

    // find current crime rate data and set it as the elevation on the map
    data.forEach(function (e) {
        var found = temparray.find(d => d.name === e.properties.name);
        if (typeof found === 'undefined') {
        } else {
            e.properties.level = found.crimerate;
        }
    });

    renderlayer();
}


/**
 * Render the GeoJsonLayer
 */
function renderlayer() {
    console.log('rendering!')

    // set a factor if the crime rate is too large
    var factor = 100;
    if ($("#dropmenu :selected").text() === 'Totaal misdrijven') {
        factor = 5;
    }

    // add the layer
    var geojsonLayer = new deck.GeoJsonLayer({

        data,
        opacity: 0.6,
        stroked: false,
        filled: true,
        extruded: true,
        wireframe: true,
        fp64: false,
        lightSettings: LIGHT_SETTINGS,

        // set the bar chart value
        getElevation: f => f.properties.level * factor,
        getFillColor: f => colorScale(f.properties.level * factor / 4000),
        getLineColor: f => [255, 255, 255],
        updateTriggers: {
            // if showLibraries changes, recompute getFillColor for each point
            getElevation: [sliderdate, slidercategory],
            getFillColor: [slidercategory, sliderdate]
        },
        transitions: {
            // transition with a duration of 3000ms
            getElevation: 1000,
            getFillColor: 1000,
        },

        pickable: true,

        // add onHover interaction, show a tooltip
        onHover: updateTooltip


    });
    Deck.setProps({
        layers: [geojsonLayer]
    });

}


const LIGHT_SETTINGS = {
    lightsPosition: [-125, 50.5, 5000, -122.8, 48.5, 8000],
    ambientRatio: 0.2,
    diffuseRatio: 0.5,
    specularRatio: 0.3,
    lightsStrength: [1.0, 0.0, 2.0, 0.0],
    numberOfLights: 2
};

const COLOR_SCALE = [
    [65, 182, 196],
    [127, 205, 187],
    [199, 233, 180],
    [237, 248, 177],
    [255, 255, 204],
    [255, 237, 160],
    [254, 217, 118],
    [254, 178, 76],
    [253, 141, 60],
    [252, 78, 42],
    [227, 26, 28],
    [189, 0, 38],
    [128, 0, 38]
];

/**
 * return the color value in COLOR_SCALE
 * @param x
 * @returns {*}
 */
function colorScale(x) {
    const i = Math.round(x * 6);
    if (x < 0) {
        return COLOR_SCALE[i] || COLOR_SCALE[0];
    }
    return COLOR_SCALE[i] || COLOR_SCALE[COLOR_SCALE.length - 1];
}


/**
 * update the information on the tooltip
 * @param x
 * @param y
 * @param object
 */
function updateTooltip({x, y, object}) {
    const tooltip = document.getElementById('tooltip');
    if (object) {
        tooltip.style.top = `${y}px`;
        tooltip.style.left = `${x}px`;
        tooltip.innerHTML = `
            <div><b>Gemeente: &nbsp;</b></div>    
            <div><div>${object.properties.name} </div></div> 
            <div><b>aantal misdrijven: &nbsp;</b></div>
            <div><div>${object.properties.level} </div></div>
        `;
    } else {
        tooltip.innerHTML = '';
    }
}